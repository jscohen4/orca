import json
def modify(file,obj,value):
	with open(file, "r") as jsonFile:
	    data = json.load(jsonFile)
	tmp = data[obj]
	data[obj] = value
	with open(file, "w+") as jsonFile:
	    json.dump(data, jsonFile, indent = 4)
	