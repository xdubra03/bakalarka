#!/usr/bin/python
import pandas as pd
file_name = "train_testa_data_new1.csv"
file_name_output = "train_testa_data_new2.csv"

df = pd.read_csv(file_name, sep=",")

# Notes:
# - the `subset=None` means that every column is used
#    to determine if two rows are different; to change that specify
#    the columns as an array
# - the `inplace=True` means that the data structure is changed and
#   the duplicate rows are gone
df.drop_duplicates(subset=None, inplace=True)

# Write the results to a different file
df.to_csv(file_name_output)