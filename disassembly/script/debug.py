import json
import itertools

f = './result/ler_ler/viewer/data.json'
f = open(f,'r')
dataset = json.load(f)

def near_origin_unit(unit):
    if 'E' in unit:
        return unit['E'][1] > 1700 and unit['E'][1] < 1760 
    else:
        return False

def touch_last(unit):
    if 'E' in unit:
        return unit['E'][1] > 1760 and unit['E'][1] < 1800
    else:
        return False


def near_origin(read):
    return any(map(near_origin_unit, read['units'])) and any(map(touch_last, read['units']))

def encode(r):
    output = r["name"];
    for unit in r['units']:
        if 'E' in unit:
            output += " U{}".format(unit['E'][1])
        else:
            output += " G{}".format(unit['G'])
    return output

print(len(dataset['reads']))

if __name__ == "__main__":
    reads = list(filter(near_origin, dataset['reads']))
    print(len(reads))
    for read in reads:
        print(encode(read))

# for r in 
# for cluster in clusters:
#     print(cluster['id'])
#     for men in cluster['members']:
#         if men['cr'].get('CP'):
#             men = men['cr']['CP']
#             print(men['contig1'])
#             print(men['contig2'])
#         elif men['cr'].get('CR'):
#             men = men['cr']['CR']
#             print(men['pos'])



