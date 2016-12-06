"""
    Script to generate vms file with non-rotating disks at turbine locations 
"""
import numpy as np
import windfarm as wf

outputfile = "farm_disks.vms"
modelfile  = "Disco_antiwarping-10mmx0.6mm.stl"

# Expand this to read from windfarm file...
farm = wf.Windfarm(Nrows=12, Ncols=6)
z = 0.1
sc = 0.010

with open(outputfile,'w') as file:
    
    # Write the header
    file.write('<?xml version="1.0" encoding="ISO-8859-1" standalone="yes"?>\n')
    file.write('\n')
    file.write('<ModelScene>\n')
    file.write('\n')
    file.write('<File>'+modelfile+'</File>\n')
    file.write('\n')
    file.write('<!-- ts 0 -->')

    file.write('<Timestep>\n')
    file.write('\n')
    cnt = 0
    for turbine in farm.turbines:

        print('Writing turbine ', cnt)
        cnt+=1

        file.write('\n')
        file.write('<Transform>\n')
        file.write('<translate>'+"{:10.4f}".format(turbine.x)+"{:10.4f}".format(turbine.y)+"{:10.4f}".format(z)+'</translate>\n')
        file.write('<rotate> 0 1 0 90 </rotate>\n')
        file.write('<scale>'+"{:10.4f}".format(sc)+"{:10.4f}".format(sc)+"{:10.4f}".format(sc)+'</scale>\n')
        file.write('<scale>0.5 0.5 0.5</scale>\n')
        file.write('</Transform>\n')
        file.write('\n')

    file.write('\n')
    file.write('</Timestep>\n')
    file.write('\n')
    file.write('</ModelScene>\n')
