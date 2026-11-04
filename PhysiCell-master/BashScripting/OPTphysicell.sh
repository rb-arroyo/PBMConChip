#!/bin/bash
#
# Runs N replicates of a PhysiCell simulation for a given set of parameters
# and copies the outputs into Results/OPT<i>.
#
# Usage: ./OPTphysicell.sh <b> <sloc> <Tpers> <Nrealizations>
#
# This script is meant to live in <PhysiCell-master>/BashScripting/.
# The PhysiCell directory is detected automatically from the script location.
# To use a different PhysiCell directory, set the PHYSICELL_DIR environment
# variable before calling the script:
#   PHYSICELL_DIR=/path/to/PhysiCell ./OPTphysicell.sh 0.15 4.0 3.0 5

echo "b: $1"
echo "sloc: $2"
echo "Tpers: $3"
echo "Nrealizations: $4"

# Set the working directory to the PhysiCell directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHYSICELL_DIR="${PHYSICELL_DIR:-$(dirname "$SCRIPT_DIR")}"
cd "$PHYSICELL_DIR" || { echo "Cannot access $PHYSICELL_DIR"; exit 1; }

make data-cleanup

rm -r "$PHYSICELL_DIR"/Results/*
rm "$PHYSICELL_DIR/BashScripting/Stop.txt"

# Replace the parameter placeholders in the XML templates
# (SESGO = bias b, VELOCIDAD = speed, TPERS = persistence time, SEMILLA = random seed)
sed 's/SESGO/'$1'/' ./config/template5.xml > ./config/template4.xml
sed 's/VELOCIDAD/'$2'/' ./config/template4.xml > ./config/template3.xml
sed 's/TPERS/'$3'/' ./config/template3.xml > ./config/template2.xml

for i in $(LANG=en_US seq 1 1 $4);
	do
	sed 's/SEMILLA/'$i'/' ./config/template2.xml > ./config/template.xml
	./project ./config/template.xml
	mkdir "$PHYSICELL_DIR/Results/OPT$i"
	cp -r "$PHYSICELL_DIR/output" "$PHYSICELL_DIR/Results/OPT$i"
done

echo "Replicates finished in PhysiCell" > "$PHYSICELL_DIR/BashScripting/Stop.txt"
