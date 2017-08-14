#!/usr/bin/env python
#-*- coding:utf-8 -*-
################################################
#File Name: download_synapse.py
#Author: C.J. Liu
#Mail: chunjie.sam.liu@gmail.com
#Created Time: Sun 30 Jul 2017 12:54:32 PM CDT
################################################

# download
import synapseclient

syn = synapseclient.Synapse()
syn.login('chunjie.sam.liu','u201012670')

# Obtain a pointer and download the data
syn7824274 = syn.get('syn7824274')

# Get the path to the local copy of the data file
filepath = syn7824274.path

# download mRNA

syn4976369 = syn.get('syn4976369')

# Get the path to the local copy of the data file
filepath = syn4976369.path
