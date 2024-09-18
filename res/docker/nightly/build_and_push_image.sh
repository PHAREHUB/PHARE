#!/usr/bin/env bash
BASEDIR=$(dirname "$0")
RELEASE=${1:-40}
REGISTRY=${2:-"129.104.6.165:32219"}
IMAGE="phare/teamcity-fedora_data"
FULL_NAME="${REGISTRY}/${IMAGE}:${RELEASE}"
$BASEDIR/build_image.sh $RELEASE $REGISTRY
docker push $FULL_NAME
