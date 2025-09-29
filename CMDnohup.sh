m=$(cat ./cmds|wc -l)
n=1
while true
do 
running=$(ps -ef|grep /usr/lib/R/bin/exec/R |grep -v grep|wc -l)
if [ $1 -gt $running ];then
eval $(cat ./cmds |head -n ${n}|tail -n 1)
n=`expr ${n} + 1`
sleep 3
else
sleep 10
fi
if [ $n -gt $m ];then
break
fi
done
