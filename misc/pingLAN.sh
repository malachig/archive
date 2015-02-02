#!/bin/bash
x=1
while [ $x -le 254 ]
do
  ping -W 1 -o -c 1 192.168.1.$x | grep icmp_seq
  x=$(( $x + 1 ))
done



