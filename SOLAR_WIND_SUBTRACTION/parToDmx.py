from parReader import Pulsar
import sys
import numpy as np

if len(sys.argv) == 1:
    print("No input .par file provided.\n")
    sys.exit()

readFile = sys.argv[1]
if readFile[-4:] != ".par":
    print("Input file is not a .par file.\n")
    sys.exit()

writeFile = ""
if len(sys.argv) == 2:
    writeFile = readFile[:-4] + ".dmx"
else:
    writeFile = sys.argv[2]
    if writeFile[-4:] != ".dmx":
        writeFile += ".dmx"

pulsar = Pulsar(readFile)

output = open(writeFile, 'w')

output.write("# " + pulsar.PSR + " dispersion measure variation\n")
output.write("# Mean DMX value = " + "{:e}".format(np.average(pulsar.DMX_ARR)) + "\n")
output.write("# Uncertainty in average DM = ########\n")
output.write("# Columns: DMXEP DMX_value DMX_var_err DMXR1 DMXR2 DMXF1 DMXF2 DMX_bin\n")

for i in range(len(pulsar.DMX_ARR)):
    output.write("{:.4f}".format(pulsar.DMXEP_ARR[i]) + " " +
                "{:e}".format(pulsar.DMX_ARR[i]) + " " +
                "{:.3e}".format(pulsar.DMX_ERROR_ARR[i]) + " " +
                "{:.4f}".format(pulsar.DMXR1_ARR[i]) + " " +
                "{:.4f}".format(pulsar.DMXR2_ARR[i]) + "  " +
                "{:.2f}".format(pulsar.DMXF1_ARR[i]) + " " +
                "{:.2f}".format(pulsar.DMXF2_ARR[i]) + " " +
                "DX" + str(i + 1).zfill(3) + "\n")

output.close()
        

