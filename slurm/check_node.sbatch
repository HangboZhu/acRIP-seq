#!/usr/bin/env bash
#SBATCH --job-name=checknode
#SBATCH --partition=core56,corexd192,c128m1T,vmcore128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/check_node.%j.out
#SBATCH --error=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/check_node.%j.err
hostname
echo "PATH=$PATH"
ls -l /usr/sbin/unsquashfs /usr/bin/unsquashfs 2>&1
command -v unsquashfs || echo NO_UNSQUASHFS
command -v mksquashfs || echo NO_MKSQUASHFS
