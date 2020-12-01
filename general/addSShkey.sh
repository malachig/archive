#!/bin/bash

E_BADARGS=65
if [ ! -n "$1" ]
then
  echo "Usage: `basename $0` malachig newssh.bcgsc.ca [username on remote host] [remote host]"
  exit $E_BADARGS
fi  

USER=$1
RHOST=$2

echo
echo ssh-keygen -t rsa [hit enter when prompted for passphrase]
if [ -e ~/.ssh/id_rsa ]
then
  echo
  echo "~/.ssh/id_ras already exists"
else
  echo
  ssh-keygen -t rsa
fi

echo
echo "ssh $USER@$RHOST mkdir -p .ssh"
ssh $USER@$RHOST mkdir -p .ssh

echo
echo "cat ~/.ssh/id_rsa.pub | ssh $USER@$RHOST 'cat >> .ssh/authorized_keys'"
cat ~/.ssh/id_rsa.pub | ssh $USER@$RHOST 'cat >> .ssh/authorized_keys'

