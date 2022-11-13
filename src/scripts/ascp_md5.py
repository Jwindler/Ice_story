#!/usr/bin/env python
# encoding: utf-8

"""
@author: ice
@contact: storyice@163.com
@file: ascp_md5.py
@time: 2022/3/1 下午5:37
@function: 根据ENA数据的CSV文件导出MD5文件和ASCP 下载文件
"""


class GenerateAscpMd5(object):
    """
    根据ENA信息，生成ascp下载文件和md5检查文件

    参数：
        ascp_path：      ENA信息文件绝对路径
        save_path：      输出文件存放的绝对路径
        openssh_path：   ASCP下载需要的openssh文件绝对路径
    """

    def __init__(self, ascp_path, save_path, openssh_path):
        """初始化类参数

        输入包含3个参数

        ascp_path：      ENA信息文件绝对路径
        save_path：      输出文件存放的绝对路径
        openssh_path：   ASCP下载需要的openssh文件绝对路径
        """
        self.ascp_path = ascp_path
        self.save_path = save_path
        self.openssh_path = openssh_path

    def generate_ascp(self):
        md5 = {}
        url = []

        with open(self.ascp_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                temp_line = line.strip().split("\t")
                one_two = len(temp_line[6])

                # 跳过开头
                if line.startswith("study_accession"):
                    continue
                else:
                    # 判断单双端数据
                    # 单端
                    if one_two == 32:
                        md5[temp_line[7].split(";")[0].split("/")[-1]] = temp_line[6]
                        url.append(temp_line[8].split(";")[0])

                    # 双端
                    else:
                        url.append(temp_line[8].split(";")[0])
                        url.append(temp_line[8].split(";")[1])
                        md5[temp_line[7].split(";")[0].split("/")[-1]] = temp_line[
                            6
                        ].split(";")[0]
                        md5[temp_line[7].split(";")[1].split("/")[-1]] = temp_line[
                            6
                        ].split(";")[1]

        # 写入md5文件
        file_temp = self.save_path + "/md5.txt"
        with open(file_temp, "x") as f:
            for key, values in md5.items():
                f.write(values + " " + key + "\n")

        # 写如ASCP文件
        file_temp = self.save_path + "/ascp.sh"
        with open(file_temp, "x") as f:
            for i in url:
                f.write(
                    "ascp -QT -l 300m -P33001 -i "
                    + self.openssh_path
                    + " era-fasp@"
                    + i
                    + " ."
                    + "\n"
                )
