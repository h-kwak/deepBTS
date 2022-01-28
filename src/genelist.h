#ifndef _GENELIST_H
#define _GENELIST_H

#include "datafile.h"
#include <vector>
#include <string>
#include <iostream>

struct genelist
{
    int n;
	std::vector<std::string> geneID;
    std::vector<std::string> geneName;
    std::vector<int> txStart;
    std::vector<int> txEnd;
    std::vector<int> cdsStart;
    std::vector<int> cdsEnd;
    std::vector<std::vector<int> > exonStart;
    std::vector<std::vector<int> > exonEnd;
    std::vector<int> exonCount;
    std::vector<std::string> chrName;
    std::vector<std::string> strand;
   	std::vector<int> exonLength;

	void parse_exon_string(std::string s,int c,std::vector<int> &v)
	{
		std::string str;
		int spacer;
		for(int i=0;i<c;++i)
		{
			spacer=s.find_first_of(',');
			str=s.substr(0,spacer);
			v.push_back(atoi(str.c_str()));
			s=s.substr(spacer+1);
		}
	};

	void load(char *filename, std::string listType="ucsc_table")
	{
		datafile df;
		if(listType!="ucsc_table") df.load(filename,0);
		else df.load(filename);
		
		if(listType=="ucsc_table")
		{
			//name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	name2
			df.get_data(geneID,0);
			df.get_data(geneName,10);
			df.get_data(chrName,1);
			df.get_data(strand,2);
			df.get_data(exonCount,7);
			df.get_data(txStart,3);
			df.get_data(txEnd,4);
			df.get_data(cdsStart,5);
			df.get_data(cdsEnd,6);
		}
		else if(listType=="bed12")
		{
			//chrom	txStart	txEnd	name	score	strand	cdsStart	cdsEnd	itemRGB	exonCount	exonSizes	exonStarts
			df.get_data(geneName,3);
			df.get_data(chrName,0);
			df.get_data(strand,5);
			df.get_data(exonCount,9);
			df.get_data(txStart,1);
			df.get_data(txEnd,2);
			df.get_data(cdsStart,6);
			df.get_data(cdsEnd,7);
			geneID=geneName;
		}
		else if(listType=="bed6")
		{
			//chrom	txStart	txEnd	name	score	strand	cdsStart	cdsEnd	itemRGB	exonCount	exonSizes	exonStarts
			df.get_data(geneName,3);
			df.get_data(chrName,0);
			df.get_data(strand,5);
			df.get_data(txStart,1);
			df.get_data(txEnd,2);
			geneID=geneName;
			exonCount.resize(geneName.size(),0);	
		}
		else if(listType=="bed3")
		{
			//chrom	txStart	txEnd
			df.get_data(chrName,0);
			df.get_data(txStart,1);
			df.get_data(txEnd,2);
			strand.resize(chrName.size(),"+");
		}

		n=chrName.size();
		
		exonStart.resize(n);
		exonEnd.resize(n);
		exonLength.resize(n);

		if(listType=="ucsc_table")
		{
			std::vector<std::string> exonStartString,exonEndString;
			df.get_data(exonStartString,8);
			df.get_data(exonEndString,9);
			
			for(int i=0;i<n;++i)
			{
				parse_exon_string(exonStartString[i],exonCount[i],exonStart[i]);
				parse_exon_string(exonEndString[i],exonCount[i],exonEnd[i]);
				for(int j=0;j<exonCount[i];++j) exonLength[i]+=exonEnd[i][j]-exonStart[i][j];
			}
		}
		else if(listType=="bed12")
		{
			std::vector<std::string> exonSizeString,exonStartString;
			df.get_data(exonSizeString,10);
			df.get_data(exonStartString,11);
				
			for(int i=0;i<n;++i)
			{
				std::vector<int> exonSize;
				std::vector<int> exonRelStart;
				parse_exon_string(exonSizeString[i],exonCount[i],exonSize);
				parse_exon_string(exonStartString[i],exonCount[i],exonRelStart);
				for(int j=0;j<exonCount[i];++j)
				{
					exonStart[i].push_back(txStart[i]+exonRelStart[j]);
					exonEnd[i].push_back(exonStart[i][j]+exonSize[j]);
					exonLength[i]+=exonSize[j];
				}
			}
		}
	};
};
#endif
