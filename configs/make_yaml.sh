
ls /home/yuansh/Desktop/Normal_stomach/ | grep _dge.csv | sed 's/_dge.csv//' | while read id
do 
cat BscModel_V4_configs.yaml | sed "s/GSE131907/${id}/g" > ${id}.yaml
done



ls /home/yuansh/Desktop/Normal_stomach/ | grep yaml | while read id 
do
python BscModel_V4.py --config /home/yuansh/Desktop/Normal_stomach/$id
done


ls /media/yuansh/14THHD/BscModel-V4/GSE148673 | grep ATC | grep yaml | while read id 
do
python BscModel_V4.py --config /media/yuansh/14THHD/BscModel-V4/GSE148673/$id
done



python BscModel_V4.py --config /home/yuansh/Desktop/Normal_stomach/Fetal-Stomach2.yaml

python BscModel_V4.py --config /home/yuansh/Desktop/GSE131907/BscModel_V4_configs.yaml

python BscModel_V4.py --config /home/yuansh/Desktop/PMID33958794/BscModel_V4_configs.yaml

python BscModel_V4.py --config /home/yuansh/Desktop/PMID33958794/BscModel_V4_configs-1.yaml
