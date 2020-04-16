import json
import itertools

f = './result/pacbio_no/viewer/data.json'
f = open(f,'r')
dataset = json.load(f)


clusters = [x for x in dataset['clusters']]


for cluster in clusters:
    print(cluster['id'])
    for men in cluster['members']:
        if men['cr'].get('CP'):
            men = men['cr']['CP']
            print(men['contig1'])
            print(men['contig2'])
        elif men['cr'].get('CR'):
            men = men['cr']['CR']
            print(men['pos'])



