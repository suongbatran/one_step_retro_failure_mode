import os
import json
import gzip
import logging
import pandas as pd
from typing import List, Dict, Any, Tuple


def setup_logging(log_file=None, log_format='%(asctime)s - %(levelname)s - %(message)s'):
    """
    Setup logging configuration.
    
    Args:
        log_file: Path to log file. If None, logs only to console.
        log_format: Format string for log messages.
    """
    if log_file:
        # Create log directory if it doesn't exist
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        logging.basicConfig(
            level=logging.INFO,
            format=log_format,
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()  # Also log to console
            ]
        )
    else:
        logging.basicConfig(
            level=logging.INFO,
            format=log_format
        )

def save_reactions_from_list(
    reactions: List[Dict[str, Any]],
    template_set: str,
    reaction_file_for_preprocessing: str,
    reaction_file_for_db: str
) -> None:
    logging.info(f"Saving unique reactions for preprocessing to "
                 f"{reaction_file_for_preprocessing}")

    with open(reaction_file_for_preprocessing, "w") as of:
        of.write("rxn_id,rxn_smiles\n")

        for reaction in reactions:
            of.write(f"{reaction['id']},{reaction['rxn_smiles']}\n")

            reaction["reaction_id"] = reaction.pop("id")
            reaction["reaction_smiles"] = reaction.pop("rxn_smiles")
            if "_id" not in reaction:
                reaction["_id"] = f"{template_set}_{reaction['reaction_id']}"
            if "template_set" not in reaction:
                reaction["template_set"] = template_set

    logging.info(f"Saving reactions for db seeding to {reaction_file_for_db}")
    with gzip.open(reaction_file_for_db, "wt", encoding="UTF-8") as zipfile:
        json.dump(reactions, zipfile)

    del reactions


def save_templates_from_list(templates: List[Dict[str, Any]], template_file: str):
    with open(template_file, "w") as of:
        for template in templates:
            if "rxn" in template:
                del template["rxn"]
            of.write(f"{json.dumps(template)}\n")


def save_templates_from_dict(
    templates: Dict[str, Dict[str, Any]],
    template_file: str,
    template_file_for_db: str
) -> None:
    templates_as_list = []

    logging.info(f"Saving templates to {template_file}")
    with open(template_file, "w") as of:
        for canon_templ, metadata in templates.items():
            assert metadata["reaction_smarts"] == canon_templ
            if "rxn" in metadata:
                del metadata["rxn"]
            of.write(f"{json.dumps(metadata)}\n")

            templates_as_list.append(metadata)

    logging.info(f"Saving templates for db seeding to {template_file_for_db}")
    with gzip.open(template_file_for_db, "wt", encoding="UTF-8") as zipfile:
        json.dump(templates_as_list, zipfile)


def load_templates_as_list(template_file: str
                           ) -> Tuple[List[Dict[str, Any]], pd.DataFrame]:
    templates = []
    template_attributes = []

    with open(template_file, "r") as f:
        for line in f:
            template = json.loads(line.strip())
            del template["references"]

            templates.append(template)
            template_attributes.append(template.get("attributes", {}))
    template_attributes = pd.DataFrame(template_attributes)

    return templates, template_attributes


def load_templates_as_dict(template_file: str) -> Dict[str, Dict[str, Any]]:
    templates = {}

    with open(template_file, "r") as f:
        for line in f:
            template = json.loads(line.strip())
            templates[template["reaction_smarts"]] = template

    return templates
