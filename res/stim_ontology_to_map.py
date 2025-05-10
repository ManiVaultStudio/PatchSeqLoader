import json
import csv

# Load JSON from file
with open('stimulus_ontology.json', 'r') as f:
    data = json.load(f)

# Open CSV for writing
with open('output.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['code', 'name'])  # Header row

    for item in data:
        code_list = None
        name = None

        for pair in item:
            if pair[0] == 'code':
                code_list = pair[1:]
            elif pair[0] == 'name':
                name = pair[1]

        # Write one row per code
        if code_list and name is not None:
            for code in code_list:
                writer.writerow([code, name])
