set -ex

echo "sonar.host.url=https://sonarcloud.io
sonar.organization=jeandet-github
sonar.projectKey=lpp_phare
sonar.projectName=Phare
sonar.projectVersion=1.0
sonar.branch.name=%branch% 
sonar.exclusions=**/subprojects/**/*,**/googletest/**/*,**/html/**/*
sonar.cfamily.build-wrapper-output=.
sonar.cfamily.threads=4
sonar.sources=.
sonar.language=cpp
sonar.cfamily.cppunit.reportsPath=.
sonar.cxx.coverage.reportPath=./coverage/coverage.xml
sonar.cfamily.gcov.reportsPath=.">%system.teamcity.build.checkoutDir%/sonar-project.properties
