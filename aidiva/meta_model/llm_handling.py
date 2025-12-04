import argparse
import gzip
import json
import logging
import os
import pandas as pd
import random
import re
import time


logger = logging.getLogger(__name__)


# generate mapping HPO -> translation
def generate_hpo2name_dict(hpo_ontology):
    hpo_mapping = dict()
    token = False

    with open(hpo_ontology, "r") as ontology:
        for line in ontology:
            if line.startswith("[Term]"):
                if token:
                    hpo_mapping[hpo_term] = dict()
                    hpo_mapping[hpo_term]["name"] = name
                    hpo_mapping[hpo_term]["synonyms"] = synonyms

                token = True
                synonyms = []
                name = ""

            elif line.startswith("id:"):
                hpo_term = line.strip().split("id: ")[1]

            elif line.startswith("name:"):
                name = line.strip().split("name: ")[1]

            elif line.startswith('synonym:'):
                synonyms.append(line.strip().split('"')[1])

            elif line.startswith("[Typedef]"):
                break

            else:
                continue

        # add information of the last entry of the file to the HPO edges file
        if token:
            hpo_mapping[hpo_term] = dict()
            hpo_mapping[hpo_term]["name"] = name
            hpo_mapping[hpo_term]["synonyms"] = synonyms

    logger.info(f"HPO to name mapping successfully loaded and prepared")
    return hpo_mapping


def get_resource_file(resource_path):
    if resource_path.startswith(".."):
        resource_file = os.path.dirname(__file__) + "/" + resource_path

    else:
        resource_file = resource_path

    return resource_file


def get_hpo_translations(hpo_terms, hpo2name_mapping):
    hpo_translations = []

    for term in hpo_terms.split(","):
        if term in hpo2name_mapping.keys():
            hpo_translations.append(hpo2name_mapping[term]["name"])

        else:
            logger.warning(f"No HPO translation found! Skip HPO term ({term})!")

    return hpo_translations


def create_llm_prompt(sex, age, hpo_terms, top_ranking_genes, model_info, internal_parameter_dict, CONSTANT_DICTIONARY, with_variant_type=True, with_genotype=True):
    # get constants
    VARIANT_CONSEQUENCES = CONSTANT_DICTIONARY["VARIANT_CONSEQUENCES"]
    CONSEQUENCE_MAPPING = CONSTANT_DICTIONARY["CONSEQUENCE_MAPPING"]
    SUPPORTED_VARIANT_TYPES = CONSTANT_DICTIONARY["SUPPORTED_VARIANT_TYPES"]
    RANDOM_SEED = CONSTANT_DICTIONARY["RANDOM_SEED"]

    hpo_ontology_f = get_resource_file(internal_parameter_dict["hpo-ontology"])
    hpo2name_mapping = generate_hpo2name_dict(hpo_ontology_f)
    hpo_translations = get_hpo_translations(hpo_terms, hpo2name_mapping)
    phenotype_information = ", ".join(hpo_translations)
    gene_list = top_ranking_genes
    list_of_gene_dicts = []

    if model_info == "rf":
        splitted_gene_list = gene_list.split(";")

    elif model_info == "eb_dom" or model_info == "eb_rec":
        splitted_gene_list = gene_list.split(";")

    for entry in splitted_gene_list:
        if model_info == "rf":
            gene_name = entry.split("(")[0]
            gene_info = entry.split("(")[1].replace(")", "")
            consequence = gene_info.split(",")[0].split(": ")[1]
            variant_type = gene_info.split(",")[1].split(": ")[1]
            rank = gene_info.split(",")[2].split(": ")[1]
            score = gene_info.split(",")[3].split(": ")[1]
            genotype = gene_info.split(",")[4].split(": ")[1]

        elif model_info == "eb_dom" or model_info == "eb_rec":
            gene_name = entry.split("(")[0]
            gene_info = entry.split("(")[1].replace(")", "")
            consequence = gene_info.split(",")[0].split(": ")[1]
            variant_type = gene_info.split(",")[1].split(": ")[1]
            rank = gene_info.split(",")[2].split(": ")[1]
            score = gene_info.split(",")[3].split(": ")[1]
            genotype = gene_info.split(",")[4].split(": ")[1]

            if "/" in gene_name:
                if len(gene_consequence.split("/")) == 2:
                    if any(consequence in gene_consequence.split("/")[0] for consequence in CONSEQUENCE_MAPPING.keys()) and not any(consequence in gene_consequence.split("/")[1] for consequence in CONSEQUENCE_MAPPING.keys()):
                        gene_info = gene_info.split("/")[0]
                        variant_type = variant_type.split("/")[0]

                    elif any(consequence in gene_info.split("/")[1] for consequence in CONSEQUENCE_MAPPING.keys()) and not any(consequence in gene_info.split("/")[0] for consequence in CONSEQUENCE_MAPPING.keys()):
                        gene_name = gene_name.split("/")[1]
                        variant_type = variant_type.split("/")[1]

                    else:
                        logger.error("No supported variant type present!!!")

            else:
                gene_name = gene_name.split("/")[0]
                variant_type = variant_type.split("/")[0]

            if "+" in variant_type:
                consequences = [min([VARIANT_CONSEQUENCES.get(consequence) if consequence in VARIANT_CONSEQUENCES.keys() else VARIANT_CONSEQUENCES.get("unknown") for consequence in variant_type.split("+")])]
                target_index = min(enumerate(consequences), key=itemgetter(1))[0]
                variant_type = variant_type.split("+")[target_index]

                if variant_type in CONSEQUENCE_MAPPING.keys():
                    variant_type = CONSEQUENCE_MAPPING[variant_type]

                else:
                    variant_type = "unspecified variant"

            else:
                if variant_type in CONSEQUENCE_MAPPING.keys():
                    variant_type = CONSEQUENCE_MAPPING[variant_type]

                else:
                    if not variant_type in SUPPORTED_VARIANT_TYPES:
                        variant_type = "unspecified variant"

        if genotype == "het":
            genotype = "heterozygous"

        elif genotype == "hom":
            genotype = "homozygous"

        list_of_gene_dicts.append({"gene": gene_name, "variant_type": variant_type, "genotype": genotype})

    random.seed(RANDOM_SEED)
    random.shuffle(list_of_gene_dicts)

    if sex != "" and str(age) != "nan":
        llm_prompt = f"A {sex} rare disease patient of age {int(age)} has the following symptoms: {phenotype_information}. "

    elif sex != "" and str(age) == "nan":
        llm_prompt = f"A {sex} rare disease patient has the following symptoms: {phenotype_information}. "

    elif sex == "" and str(age) != "nan":
        llm_prompt = f"A rare disease patient of age {int(age)} has the following symptoms: {phenotype_information}. "

    else:
        llm_prompt = f"A rare disease patient has the following symptoms: {phenotype_information}. "

    llm_prompt += "A causal variant in which of the following candidate genes would best explain these symptoms? "

    if with_variant_type:
        if with_genotype:
            causal_genes = []

            for entry in list_of_gene_dicts:
                if entry['genotype'] == "homozygous":
                    causal_genes.append(str(f"{entry['gene']} ({entry['genotype']} {entry['variant_type']})"))

                else:
                    causal_genes.append(str(f"{entry['gene']} ({entry['variant_type']})"))

            llm_prompt += "For each candidate gene the type of the variant is given in brackets. Furthermore the genotype of the variant is included in the brackets if it is a homozygous variant. "
            llm_prompt += f"Candidate genes: {', '.join(causal_genes)}"

        else:
            causal_genes = []

            for entry in list_of_gene_dicts:
                causal_genes.append(str(f"{entry['gene']} ({entry['variant_type']})"))

            llm_prompt += "For each candidate gene the type of the variant is given in brackets. "
            llm_prompt += f"Candidate genes: {', '.join(causal_genes)}"

    else:
        causal_genes = []

        for entry in list_of_gene_dicts:
            causal_genes.append(str(f"{entry['gene']}"))

        llm_prompt += f"Candidate genes: {', '.join(causal_genes)}"

    return llm_prompt


def call_llm_api(client, prompt, llm_instructions, model_id, llm_api, use_random_seed=False):
    # Define your messages for the chat
    messages = [
                {"role": "system", "content": llm_instructions},
                {"role": "user", "content": prompt}
               ]

    if llm_api == "OPENAI" or llm_api == "LOCAL":
        # Make a request to the OpenAI API using the chat endpoint
        response = client.chat.completions.create(model=model_id, messages=messages, temperature=0.1)

    #elif llm_api == "MISTRALAI":
    #    response = client.chat.complete(model=model_id, messages=messages, temperature=0.1)

    else:
        raise SystemExit(f"Unsupported API type ({llm_api})!!!")

    num_retries = 3
    for i in range(num_retries):
        try:
            # Extract the text of the response
            answer = response.choices[0].message.content.strip()
            finish_reason = response.choices[0].finish_reason
            used_tokens = response.usage

            # Print the response
            logger.debug("Response from OpenAI API:")
            logger.debug(answer, "\n\n")
            logger.debug("Tokens used:", used_tokens)
            logger.debug("Finish reason:", finish_reason, "\n\n")

            if "```json" in answer:
                answer = answer.replace("```json", "```")

            cleaned_answer = answer.split("```")
            logger.debug(answer, "\n\n")
            logger.debug("Cleaned answer:\n", cleaned_answer)

            if len(cleaned_answer) > 1:
                json_answer = cleaned_answer[1]

            else:
                json_answer = cleaned_answer[0]

            logger.debug("JSON answer:\n", json_answer)

            if json_answer.startswith("{"):
                if not json_answer.endswith("}"):
                    json_answer = json_answer + "}"

            # find a fix for the problem that sometimes before and/or after the json data additional sentences can occur depending on the used LLM
            list_of_answers = json.loads(json_answer)

        except Exception as e:
            if i < num_retries - 1: # i is zero indexed
                logger.warning(f"Caught {type(e)}: {e}")
                logger.warning(f"Retry in 30 seconds to get answer from LLM! ({i + 1}  out of {num_retries} tries)")
                time.sleep(30)
                continue

            else:
                logger.error(f"Failed {num_retries} times! Execution stopped please investigate the problem manually!")
                raise

        break

    if not isinstance(list_of_answers, list):
        if isinstance(list_of_answers, dict):
            list_of_answers = [list_of_answers]

    results = pd.DataFrame(list_of_answers)

    return results


def check_gene_rank(row):
    causal_genes = row["causal genes"]

    if not row.empty:
        llm_gene_1 = str(llm_response["1st ranked Gene"].values[0]).upper()
        llm_gene_2 = str(llm_response["2nd ranked Gene"].values[0]).upper()
        llm_gene_3 = str(llm_response["3rd ranked Gene"].values[0]).upper()

    else:
        logger.warning("Empty result list!")
        llm_gene_1 = "nan"
        llm_gene_2 = "nan"
        llm_gene_3 = "nan"

    if llm_gene_1 in causal_genes:
        llm_causal_rank = "1"

    elif llm_gene_2 in causal_genes:
        llm_causal_rank = "2"

    elif llm_gene_3 in causal_genes:
        llm_causal_rank = "3"

    else:
        llm_causal_rank = "-1"

    return llm_causal_rank


# for debugging
def main(infile, outfile, result_path, rank):
    pass


## The possibility to directly run the script is mainly meant for testing und debugging
if __name__ == "__main__":
    pass
