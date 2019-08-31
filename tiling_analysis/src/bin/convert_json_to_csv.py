import sys
import json
import csv
if __name__ == "__main__":
    ARGS = sys.argv;
    with open(ARGS[1], 'r') as file:
        input = json.load(file)
        with open(ARGS[2], 'w') as output:
            fieldnames = ["source", "target", "weight"]
            writer = csv.DictWriter(output, fieldnames=fieldnames)
            writer.writeheader()
            for edge in input["edges"]:
                writer.writerow({"source":edge[0], "target":edge[1], "weight":edge[2]})
        with open(ARGS[3], 'w') as output:
            fieldnames = ["id", "color"]
            writer = csv.DictWriter(output, fieldnames=fieldnames)
            writer.writeheader()
            for node in input["nodes"]:
                writer.writerow({"id":node[0], "color":node[1]})
