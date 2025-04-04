#!/bin/bash

for ((a=2000; a<=2000; a++))

do

        echo "#!/bin/bash"  >  ma-$a.bash
        echo "#PBS -l walltime=168:00:00" >> ma-$a.bash
        echo "#PBS -q old" >> ma-$a.bash
        echo "#PBS -l nodes=1:ppn=12" >> ma-$a.bash                             
        echo "#PBS -m bea" >> ma-$a.bash
        echo "#PBS -M luis.serranof@estudiante.uam.es" >> ma-$a.bash                                  
        echo "#PBS -e /home/alvaro/javier/codificacion_1/output-E-"$a >> ma-$a.bash        
        echo "#PBS -o /home/alvaro/javier/codificacion_1/output-O-"$a >> ma-$a.bash
        echo "/usr/local/bin/matlab -nosplash -nodesktop nodisplay -r \" cd /home/alvaro/javier/codificacion_1/;  aux_function("$a")\"" >> ma-$a.bash        

        qsub ma-$a.bash

done
