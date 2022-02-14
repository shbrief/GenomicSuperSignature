#!/bin/bash

#' To be used after a Bioconductor release
#'
#' Use `git remote -v` to check if you have the correct origin and upstream
#' remotes 
#'
#' If NOT, use
#      `git remote add upstream git@git.bioconductor.org:packages/PACKAGE`
#' If your origin is not GitHub, use
#'     `git remote origin set-url git@github.com:USER/PACKAGE`
#'
#' param RELEASE the release branch name e.g., RELEASE_3_14
#'
#' usage:
#' cd PACKAGE
#' ./update_release.sh RELEASE_3_14
#' ./update_release.sh RELEASE_3_14 main

RELEASE=$1
DEFAULT=${2:-master}

git fetch --all
git pull upstream $DEFAULT
git push origin $DEFAULT
git checkout -b $RELEASE upstream/$RELEASE
git push origin $RELEASE
git checkout $DEFAULT
