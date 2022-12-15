import pandas as pd

# Load data from excel file
path = '/Users/anders/Documents/Master/data/SOLTI-1402 CORALLEEN for UiO.xlsx'
data = pd.read_excel(path)

# Access all data in a specific column: data['Ki67 (%)']
# Subset data to only contain rows with specific column entries: data[data['Trial Arm Neo'] == 'Letro+Ribo']
# here need to filter out nan for Y values: for i in range(51): print(np.isnan(data['Ki67 (%)'][i]))

# df_Ki67 = data[np.isnan(data['Ki67 (%)']) == False]
df_Ki67 = data[data['Ki67 (%)'].notna()]
df_Ki67 = df_Ki67[df_Ki67['Trial Arm Neo'] == 'Letro+Ribo']
gene_columns = data.keys()[40:]
df_Ki67
df_Ki67_dropNAN = df_Ki67.dropna(axis='index', subset= gene_columns)

# I wrote some small code to get gene expression data for some specific genes
# for all timepoints and all patients which looks like this  (it is not very elegant)

patient_ids = list(set(data['UniqueID']))
# Be aware, that the order of the patient ids will not be the same as in the "data" df, since "set" does not preserve ordering
num_patients = len(patient_ids)

expression = {}
gene_ids = ['MYC', 'RB1', 'CCND1', 'CCNE1', 'ESR1', 'CDKN1A']
data.keys()[40:]

for gene_id in gene_ids:
    expression[gene_id] = np.full(shape=(num_patients, 3), fill_value=np.nan)
    for idx, patient_id in enumerate(patient_ids):
        patient_data = data[data['UniqueID'] == patient_id]
        t_scr = patient_data[patient_data['timepoint'] == 'SCR']
        t_sur = patient_data[patient_data['timepoint'] == 'SUR']
        t_2weeks = patient_data[patient_data['timepoint'] == '3weeks']
        if len(t_scr[gene_id]) > 0:
            expression[gene_id][idx, 0] = t_scr[gene_id].values[0]
        if len(t_2weeks[gene_id]) > 0:
            expression[gene_id][idx, 1] = t_2weeks[gene_id].values[0]
        if len(t_sur[gene_id]) > 0:
            expression[gene_id][idx, 2] = t_sur[gene_id].values[0]

