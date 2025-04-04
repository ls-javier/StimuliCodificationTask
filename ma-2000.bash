#!/bin/bash
#PBS -l walltime=168:00:00
#PBS -q old
#PBS -l nodes=1:ppn=12
#PBS -m bea
#PBS -M luis.serranof@estudiante.uam.es
#PBS -e /home/alvaro/javier/codificacion_1/output-E-2000
#PBS -o /home/alvaro/javier/codificacion_1/output-O-2000
/usr/local/bin/matlab -nosplash -nodesktop nodisplay -r " cd /home/alvaro/javier/codificacion_1/;  aux_function(2000)"
