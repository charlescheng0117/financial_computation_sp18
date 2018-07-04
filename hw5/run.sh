FILE=1.in

echo "Compiling source code"
make
echo "Done"
echo "Demo Binomial Tree methods"
time ./hw5_binomial < $FILE
echo "Demo Monte Carlo methods"
time ./hw5_monte < $FILE
echo "Finished."
make clean
