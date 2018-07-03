echo "Compiling source code"
make
echo "Done"
echo "Demo Binomial Tree methods"
time ./hw5_binomial < 1.in
echo "Demo Monte Carlo methods"
time ./hw5_monte < 1.in
echo "Finished."
make clean
