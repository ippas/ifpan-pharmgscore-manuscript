import pandas as pd
import re

source_dir = "../data/lkps/"
file_with_list_of_drugs = "../data/ukb-mdd-phenotypes.csv" 

list_of_drugs = pd.read_csv(file_with_list_of_drugs, usecols=[1])['term'].tolist()
bnf_df = pd.read_csv(source_dir + "bnf_lkp.csv")
read_disease_df = pd.read_csv(source_dir + "read_v2_lkp.csv")
read_drug_df = pd.read_csv(source_dir + "read_v2_drugs_lkp.csv")
dmd_df = pd.read_csv(source_dir + "dmd_lkp.csv", dtype={'concept_id': str, 'term': str,})
ctv3_df = pd.read_csv(source_dir + "read_ctv3_lkp.csv")

print(list_of_drugs)

def convert_bnf_code(bnf_code):
    converted_bnf_code = bnf_code[0] + bnf_code[1] + '.' + bnf_code[2] + bnf_code[3] + '.' + bnf_code[4] + bnf_code[5] + ".00.00"
    return converted_bnf_code

def find_drug(field, drugs_list):
    for drug in drugs_list:
        if drug in field:
            return drug
        if drug.lower() in field:
            return drug
        if drug.upper() in field:
            return drug
    return None

def edit_name(drugname):
    return drugname.split('_')[0].split(' ')[0].strip()

def set_term(x, drug_brand_name_list):
    if x is None:
        return ""
    for i in range(len(drug_brand_name_list)):
        if x.upper() in drug_brand_name_list[i][1]:
            return drug_brand_name_list[i][0]
        if x.lower() in drug_brand_name_list[i][1]:
            return drug_brand_name_list[i][0]
        if x in drug_brand_name_list[i][1]:
            return drug_brand_name_list[i][0]

def extract_dose(trade_name):
    # Match x mg/ x mL
    mg_ml_matches = re.findall(r'(\d+)\s?mg/(\d+)\s?ml', trade_name, re.IGNORECASE)
    if mg_ml_matches:
        mg_value, ml_value = mg_ml_matches[0]
        return f"{mg_value} mg/{ml_value} ml"
    
    # Match mg/mL
    mg_ml_matches = re.findall(r'(\d+)\s?mg/mL', trade_name, re.IGNORECASE)
    if mg_ml_matches:
        # Return the first mg/mL match (assuming there's only one)
        return f"{mg_ml_matches[0]} mg/mL"
    
    # Match mg
    mg_matches = re.findall(r'(\d+)\s?mg\b', trade_name, re.IGNORECASE)
    if mg_matches:
        # Sum all mg matched
        total_mg = sum(int(match) for match in mg_matches)
        return f"{total_mg} mg"
    
    return None
        
def extract_tablets(trade_name):
    # Match tablets
    tablet_matches = re.findall(r'(\d+)\s?(?:tabs?|tablets?|tab)\b', trade_name, re.IGNORECASE)
    if tablet_matches:
        return int(tablet_matches[0])
    return None

bnf_df['Term'] = bnf_df['BNF_Chemical_Substance'].apply(lambda x: find_drug(str(x), list_of_drugs))
bnf_df = bnf_df.dropna(subset=['Term'])
bnf_df['Converted_BNF_Code'] = bnf_df['BNF_Presentation_Code'].apply(lambda x: convert_bnf_code(str(x)))
bnf_df['Brand_Name'] = bnf_df['BNF_Presentation'].apply(lambda x: edit_name(str(x)))

bnf_df[['Converted_BNF_Code', 'Term', 'BNF_Presentation_Code', 'BNF_Presentation', 'Brand_Name']].to_csv('../data/bnf_brand.csv', index=False)
print("bnf df done")


read_disease_df['term'] = read_disease_df['term_description'].apply(lambda x: find_drug(str(x), list_of_drugs))
read_disease_df = read_disease_df.drop_duplicates(subset=['term'])

read_disease_df[['term','read_code', 'term_description']].to_csv('../data/read_name_dis.csv', index=False)
print("read disease df done")

drug_brand_name_list = list(zip(bnf_df['Term'].tolist(), bnf_df['Brand_Name'].tolist()))
drug_brand_name_list = list(set(drug_brand_name_list))
drug_brand_name_list = [pair for pair in drug_brand_name_list if pair[1] != 'Gx']

brand_names_list = bnf_df['Brand_Name'].drop_duplicates().tolist()

read_drug_df['brand_name'] = read_drug_df['term_description'].apply(lambda x: find_drug(str(x), brand_names_list))
read_drug_df['term'] = read_drug_df['brand_name'].apply(lambda x: set_term(edit_name(str(x)), drug_brand_name_list))
read_drug_df = read_drug_df.dropna(subset=['brand_name']).dropna(subset=['term'])
read_drug_df['dose'] = read_drug_df['term_description'].apply(extract_dose)
read_drug_df['quantity'] = read_drug_df['term_description'].apply(extract_tablets)

read_drug_df[['read_code','brand_name', 'term_description', 'term', 'quantity', 'dose']].to_csv('../data/read_name_drug.csv', index=False)
print("read drug df done")

list_of_keywords = ['poisoning', 'adverse', 'level', 'measurment', 'reaction', 'contraindicated', 'refused', 'urine', 'enuresis', 'overdose', 'allergy']
ctv3_df['brand_name'] = ctv3_df['term_description'].apply(lambda x: find_drug(str(x), brand_names_list)) 
ctv3_df['term'] = ctv3_df['brand_name'].apply(lambda x: set_term(edit_name(str(x)), drug_brand_name_list))
ctv3_df = ctv3_df.dropna(subset=['brand_name', 'term'])
filtered_ctv3_df = pd.DataFrame()
for index, row in ctv3_df.iterrows():
    if not any(keyword in row['term_description'].lower() for keyword in list_of_keywords):
        filtered_ctv3_df = filtered_ctv3_df.append(row)
filtered_ctv3_df['dose'] = filtered_ctv3_df['term_description'].apply(extract_dose)
filtered_ctv3_df['quantity'] = filtered_ctv3_df['term_description'].apply(extract_tablets)
print(filtered_ctv3_df)
filtered_ctv3_df[['read_code', 'brand_name', 'term_description', 'term', 'dose', 'quantity']].to_csv('../data/ctv3_drug.csv', index=False)
print("ctv3 done")

dmd_df = dmd_df.rename(columns={'concept_id': 'dmd_code', 'term': 'presentation'})
dmd_df['brand_name'] = dmd_df['presentation'].apply(lambda x: find_drug(str(x), bnf_df['Brand_Name'].tolist()))
dmd_df['term'] = dmd_df['brand_name'].apply(lambda x: set_term(edit_name(str(x)), drug_brand_name_list))
dmd_df = dmd_df.dropna(subset=['brand_name']).dropna(subset=['term'])
dmd_df['dose'] = dmd_df['presentation'].apply(extract_dose)
dmd_df['tablets'] = dmd_df['presentation'].apply(extract_tablets)

dmd_df[['dmd_code', 'brand_name', 'presentation', 'term', 'quantity', 'dose']].to_csv('../data/dmd_name.csv', index=False)
print("dmd drug df done")
