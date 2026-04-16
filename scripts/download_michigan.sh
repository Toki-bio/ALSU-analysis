#!/bin/bash
set -e
echo "Downloading 27 files from Michigan Imputation Server..."

curl -L https://imputationserver.sph.umich.edu/share/results/2a87b30430bbd44939595486534fcb9df55426fd81789165d9077e178daccdff/chr_1.zip -o chr_1.zip
echo "1/27 chr_1.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/e50af01426196de2556bdcce7a4412b1d94d62d8a5fcf06e9d81cefa6f1ccabb/chr_10.zip -o chr_10.zip
echo "2/27 chr_10.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/d7ede4954d11b2e57880f8570b5939dfc163e4109d32897f4f6e40d22bcc75f7/chr_11.zip -o chr_11.zip
echo "3/27 chr_11.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/8c20885eef6e66c5656d19663e30772714b4bba31509a7d66dc2dc81bbf4f5f7/chr_12.zip -o chr_12.zip
echo "4/27 chr_12.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/c469a0d38d7e292f25af6c69a295365fe75910c3bd7f82065f98e0ef3bf61194/chr_13.zip -o chr_13.zip
echo "5/27 chr_13.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/5b107fa2f3eb5eb1a460e65ff1bc9e7fb74b6a9d8e1e976b986bb1aba11c3757/chr_14.zip -o chr_14.zip
echo "6/27 chr_14.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/82469309f2752b1dce8bbdff97de4c373be06e3eeef5d45c2e48c9febab1b947/chr_15.zip -o chr_15.zip
echo "7/27 chr_15.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/9171ef87052810bb33ee58a1798c8fdbf9e0606a28dfef2450fc5870394b14d2/chr_16.zip -o chr_16.zip
echo "8/27 chr_16.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/8c497ee4bea8569f9b87d334c97022bab7ddaf97f86a23932b5d28f5cb20c85d/chr_17.zip -o chr_17.zip
echo "9/27 chr_17.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/4d78b2e898814c3135a75bccd005e2f33365406a9ac0a449a7e6845dc8fe8043/chr_18.zip -o chr_18.zip
echo "10/27 chr_18.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/47566eca76332d15b02cdf8eab040ee5bc968e9fe99a088a8a2090cd288366f1/chr_19.zip -o chr_19.zip
echo "11/27 chr_19.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/ec0f1ecdd1ae72b052894113a2e77f09f8143ac2d6d6e08b6d0e78c02b62bbbf/chr_2.zip -o chr_2.zip
echo "12/27 chr_2.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/e7e3974d1fb72c7d46f700582645db69da15b4471002887a7a3f941398dd7160/chr_20.zip -o chr_20.zip
echo "13/27 chr_20.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/674be72e77b0690630a9a25f53b2e57a083d93a7f2ad5874c6e42be8c8ef6727/chr_21.zip -o chr_21.zip
echo "14/27 chr_21.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/9bb42010b77064e9e0c8b2d96d6d1fbce359a91d3636addebd815fb6f344da57/chr_22.zip -o chr_22.zip
echo "15/27 chr_22.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/b33f936ad78fa65a3d1d6d7f977c5b5d2836eba61d08cee0b445ae97e47ee887/chr_3.zip -o chr_3.zip
echo "16/27 chr_3.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/9c756df94ae17767799cdb4e422f16bba0e141c44eb788c33ba3705780d71082/chr_4.zip -o chr_4.zip
echo "17/27 chr_4.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/00add947fbf01535557af01de599813357bb4bb102ffcce528175d5c03051be8/chr_5.zip -o chr_5.zip
echo "18/27 chr_5.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/56ae398622942bb46c32f45363771a820b31990d5cff703d43b360ef814f990e/chr_6.zip -o chr_6.zip
echo "19/27 chr_6.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/d567576ebc4ff599b584002f33d7e7274d41a08e97923ef2a0ef5fbccb883e0b/chr_7.zip -o chr_7.zip
echo "20/27 chr_7.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/786395b7feab0acfb8a5f1982c5646fadde2b0264acbf3298e44147a67c6d9d3/chr_8.zip -o chr_8.zip
echo "21/27 chr_8.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/b9cb5a89a9ffc65d3dadc30732f1cefb8b6b0f4cca50f38876c0636964fe5401/chr_9.zip -o chr_9.zip
echo "22/27 chr_9.zip done"

curl -L https://imputationserver.sph.umich.edu/share/results/b6f88356af6ca105e8170aa9baf694739683793ef64633db611fe7d02f25442f/qc_report.txt -o qc_report.txt
echo "23/27 qc_report.txt done"

curl -L https://imputationserver.sph.umich.edu/share/results/9ca68f9835d0b92461e0f68e5c585fd3db30567fe8106ba37bf67a65fd04979c/quality-control.html -o quality-control.html
echo "24/27 quality-control.html done"

mkdir -p statistics
curl -L https://imputationserver.sph.umich.edu/share/results/5d697668a3cd0d5bd62b728e5ab61bec415d8d630f632b570beee788df15efe8/statistics/chunks-excluded.txt -o statistics/chunks-excluded.txt
echo "25/27 chunks-excluded.txt done"

curl -L https://imputationserver.sph.umich.edu/share/results/1194b94cb7852fd808ff6258f5ef4e747c030b3cc1ae6c4a02eb799f6919efc2/statistics/snps-excluded.txt -o statistics/snps-excluded.txt
echo "26/27 snps-excluded.txt done"

curl -L https://imputationserver.sph.umich.edu/share/results/35e268282089532828d645e5f47b3f2f582c065d6d1bb5fd99011538214db8fb/statistics/snps-typed-only.txt -o statistics/snps-typed-only.txt
echo "27/27 snps-typed-only.txt done"

echo ""
echo "All 27 files downloaded."
