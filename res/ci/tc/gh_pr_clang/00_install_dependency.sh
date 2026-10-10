set -ex

ulimit -a # print system resource limits
dnf install -y cppcheck cppcheck-htmlreport  
# this is for TC errors
git config --global --add safe.directory $PWD
lscpu
