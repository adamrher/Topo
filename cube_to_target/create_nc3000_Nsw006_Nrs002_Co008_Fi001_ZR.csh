#!/bin/csh -f

setenv tag "_nc3000_Nsw006_Nrs002_Co008_Fi001_ZR"

cd /project/amp/juliob/topo_git/Topo/cube_to_target/

touch start_cube_to_t.jnk

cp cube_to_target.nl${tag} cube_to_target.nl

./cube_to_target

##mv fort.911 Paint${tag}.dat 
mv TEST.dat TEST${tag}.dat

exit
