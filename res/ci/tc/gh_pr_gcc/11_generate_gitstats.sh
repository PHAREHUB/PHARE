set -ex

mkdir stats 

# stops error "fatal: unsafe repository"
git config --global --add safe.directory $PWD

gitstats %system.teamcity.build.checkoutDir% ./stats
