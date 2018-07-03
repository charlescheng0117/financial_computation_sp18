FILE=4.in

echo "Compiling source code"
make
echo "Done"
echo "Demo Binomial Tree methods"
time ./hw5_binomial < $FILE
echo "Demo Monte Carlo methods"
time ./hw5_monte < $FILE
echo "Demo for bonus 1"
time ./hw5_bn1 < $FILE
echo "Finished."
make clean
