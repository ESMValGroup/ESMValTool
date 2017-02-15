#!/bin/bash
# Script for removing .svn from a directory tree 

# MODIFICATION HISTORY:
# 20121121-A_GoKl: written 

  WD=$(pwd)  #/tmp
  temp1=${WD}/temp1_rm_svn_$(whoami)_$$

  if [ -e ${temp1} ]; then rm ${temp1}; fi

  if [ $# -gt 1 ]; then    #arguments given
    echo 'Error: wrong number of arguments'
    echo ' process current directory tree: rm_svn'
    echo ' process $DIR: rm_svn $DIR' 
  fi
  if [ $# -eq 0 ]; then
    DIR=$(pwd)
  fi 
  if [ $# -eq 1 ]; then
    DIR=$1
  fi

#  find $(pwd) | xargs grep $1
  echo "$0 PID: $$"
  echo "Removing .svn folders from ${DIR} and subdirectories"
  find ${DIR} -type d > ${temp1}
  cat ${temp1} |
   while read line
    do        #loop over directories
      if [ -d ${line}/.svn ]; then
        echo "removing ${line}/.svn" 
        rm -rf ${line}/.svn
      fi 
    done      #end loop over directories

  cd ${WD}

  if [ -e ${temp1} ]; then rm ${temp1}; fi

