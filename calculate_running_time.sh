# Global variable

#########################################################
wkd=/media/yuansh/My\ Passport/科研项目/深度学习模型/CapsuleNet预测肿瘤细胞
#########################################################
echo running GSE131907
start=$(date +%s)

python BscModel_V4.py --config /home/yuansh/Desktop/GSE131907/BscModel_V4_configs.yaml

end=$(date +%s) 
take=$(( end - start ))
h=$(( take / 3600 ))
m=$((( take - 3600 * h ) / 60 ))
s=$((( take - 3600 * h - 60 * m )))
echo Running Completed，using: ${h} H ${m} M ${s} S

#########################################################

start=$(date +%s)

python BscModel_V4.py --config /home/yuansh/Desktop/GSE132465/BscModel_V4_configs.yaml

end=$(date +%s) 
take=$(( end - start ))
h=$(( take / 3600 ))
m=$((( take - 3600 * h ) / 60 ))
s=$((( take - 3600 * h - 60 * m )))
echo Running Completed，using: ${h} H ${m} M ${s} S

#########################################################

start=$(date +%s)

python BscModel_V4.py --config /home/yuansh/Desktop/PRJNA591860/BscModel_V4_configs.yaml

end=$(date +%s) 
take=$(( end - start ))
h=$(( take / 3600 ))
m=$((( take - 3600 * h ) / 60 ))
s=$((( take - 3600 * h - 60 * m )))
echo Running Completed，using: ${h} H ${m} M ${s} S
#########################################################
