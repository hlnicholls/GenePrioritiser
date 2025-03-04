import pickle
import sys
import pandas as pd
import os
CURRENT_DIR = os.getcwd()
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
sys.path.append(CONFIG_DIR)
import config
script_dir = config.eda_script_path
sys.path.append(script_dir)
from updated_MissForest_Algorithm import MissForest

seed = 0
model_path = config.best_model_pkl_file

# Load the model
with open(model_path, 'rb') as file:
    best_model = pickle.load(file)

training_data = pd.read_csv(config.cleaned_training_data_ml)
columns_to_remove = ['Gene', 'label_encoded']
required_columns = [col for col in training_data.columns if col not in columns_to_remove]

print('Prioritising all GWAS genes...')
dataset_unknown = pd.read_csv(config.genes_to_prioritise)
data = dataset_unknown[required_columns]

# Impute missing values
imputer = MissForest(random_state=seed)
df3 = imputer.fit_transform(data)
X2 = pd.DataFrame(df3, index=data.index, columns=data.columns)

# Make predictions
predictions = best_model.predict(X2)

# Get predicted probabilities for each class
probabilities = best_model.predict_proba(X2)

# Creating a DataFrame for probabilities
prob_df = pd.DataFrame(probabilities, columns=[f'Probability_{cls}' for cls in best_model.classes_])

# Adding genes and their predicted labels to the probabilities DataFrame
output = pd.DataFrame(dataset_unknown['Gene'], columns=['Gene'])
output['best_model_label'] = predictions

# Concatenate the gene labels and probabilities
final_output = pd.concat([output, prob_df], axis=1)

# Save to CSV
final_output.to_csv(config.best_model_predictions, index=False)


####################################

print('Prioritising all genes...')
dataset_unknown = pd.read_csv(config.all_genes_all_features_unprocessed)
dataset_unknown.set_index('Gene', inplace=True)
data = dataset_unknown[required_columns]

# Impute missing values
imputer = MissForest(random_state=seed)
df3 = imputer.fit_transform(data)
X2 = pd.DataFrame(df3, index=data.index, columns=data.columns)

# Getting data for streamlit app:
X2.to_csv(config.all_imputed_features, index=True)

# Predictions
model = best_model
predictions = model.predict(X2)

# Predicted probabilities for each class
probabilities = model.predict_proba(X2)

# Creating a DataFrame for probabilities
prob_df = pd.DataFrame(probabilities, columns=[f'Probability_{cls}' for cls in model.classes_])

# Preparing the output
dataset_unknown = pd.read_csv(config.all_genes_all_features_unprocessed)
output = pd.DataFrame(dataset_unknown['Gene'], columns=['Gene'])
output['catboost_label'] = predictions

# Concatenate the gene labels and probabilities
final_output = pd.concat([output, prob_df], axis=1)

# Save to CSV
final_output.to_csv(config.best_model_predictions_all_genes, index=False)
