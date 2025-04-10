## Releasing a new version of PHARE

- [ ] Review file: res/cmake/release.cmake
- [ ] Update to use latest (if safe) versions
- [ ] PR updates to master
- [ ] Cut branch from master as $MAJ.$MIN.x if required
- [ ] Write/commit/push release.ver.txt
        echo $MAJ.$MIN.$PATCH > release.ver.txt
        git commit -m "releasing $MAJ.$MIN.$PATCH"
        git push origin $MAJ.$MIN.x
- [ ] Release on github as $MAJ.$MIN.$PATCH
