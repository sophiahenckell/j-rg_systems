#!/bin/sh


sed -i 's/Power/pow/g' $1
sed -i 's/Cos/cos/g' $1
sed -i 's/Sin/sin/g' $1
sed -i 's/\\\[//g' $1
sed -i 's/\]//g' $1
sed -i 's/Lambda/lambda0/g' $1
sed -i 's/Theta/theta/g' $1
sed -i 's/Sqrt/sqrt/g' $1
sed -i 's/OmegaA/omegaA/g' $1
sed -i 's/OmegaB/omegaB/g' $1
