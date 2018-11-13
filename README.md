# Splice4k

**2018.10.17** - Splice4k最终算是比较成熟了，目前在功能的使用方法上也算达到一个平衡

**2018.11.13** - 鉴于老师的要求，现在转向Python重构整个流程，所以kotlin-jvm这版，弃坑啦，弃坑啦！！好处就是，可以在pysam的基础上解决几个老大难的问题了，比如不同index的支持问题

---

## 说

Splice4k 以kotlin语言编写，添加gradle和maven两种包管理框架的支持，可编译为jar包，基于java8+环境运行。

Splice4k依赖了如下框架：
- [clikt](https://github.com/ajalt/clikt) 
	- 用以处理命令行参数
	- Handling command line parameters
- [htsjdk](https://github.com/samtools/htsjdk) 
	- 用以处理BAM/SAM格式的文件
	- Handling with BAM/SAM format files
- [log4j](https://github.com/apache/log4j) 
	- 用以处理日志输出
	- Handling log output
- [progressbar](https://github.com/ctongfei/progressbar) 
	- 用以展示进度
	- Display progress

算是正经搞得第一个程序了吧
## 目的
该软件主要是
- 从RNA-seq以及SMRT-seq (PacBio、Oxford Nanopore)数据中识别可变剪接事件
- 从BAM/SAM文件中提取splice junctions
- 从SMRT-seq数据中组装转录本

主要支持几种不同的输入文件
#### 输入文件 Input file format

- STAR SJ.out.tab
- gmap align output file (gmap with -A, but without -f parameters)
- BAM/SAM
- splice junctions file extracted by Splice4k

#### 参考文件 Reference file format
仅在ensembl官方以及NCBI的文件格式上测试过。
Only tested on official ensembl/NCBI files.

- gff3
- gtf

## 安装 Installation
1. 源码 Source code
	```bash
	cd Splice4k
	// maven packaging, jar file under ./target/
	mvn package
	
	// gradle packaging, jar file under ./build/libs
	gradle build	
	```
2. jar file
	```bash
	java -jar Splice4k-x.x.x.jar
	```

## 使用 Usage
共有四种模式
4 different mode

```bash
Usage: parameters [OPTIONS] COMMAND [ARGS]...

Options:
  -v, --version  version
  -h, --help     Show this message and exit

Commands:
  extract  Extract splice junctions from BAM/SAM files
  sgs      Identify alternative splicing events from RNA-seq
  sms      Identify alternative splicing events from SMRT-seq
  iso      Construct Isoforms through SMRT-seq data
```

#### 1. Extract splice junctions
```bash
Usage: extract [OPTIONS]

  Extract splice junctions from BAM/SAM files

Options:
  -i, --input PATH   Path to input BAM/SAM file
  -o, --output PATH  Path to output file
  -c, --count INT    Filter low abundance junctions [default: 0]
  -h, --help         Show this message and exit
```

#### 2. Identify alternative splicing events based on RNA-seq data
```bash
Usage: sgs [OPTIONS] [INPUT]...

  Identify alternative splicing events from RNA-seq

Options:
  -b, --bam PATH               			Path to BAM/SAM files, or the directory
                               contains BAM/SAM files. [default: current
                               running directory] 
                               - If input files are BAM/SAM
                                 files, this parameter won't work 
                               - If specified path to directory contains                                  BAM/SAM files corresponding to STAR  									 SJ.out.tab files, this program will 								     auto match those files 
                               - If specified BAM/SAM file with this                                      parameter, then this program will                                        calculate PSI of IR using this file
  -r, --reference PATH         		    Path to reference file
  -o, --output PATH           		   Path to output file
  -c, --count INT               		Filter low abundance junctions [default: 3]
  -e INT                       			    The error to identify whether AS event exists
                               [default: 3bp]
  -p, --process INT            		      Number of processes to use
  --overlap-exon-intron FLOAT        Minimal overlap level between exon with intron
                   			required for intron retention identification 
                   			[default: 90.0]
  -v, -verbose                 			Enable detailed messages
  -h, --help                       		       Show this message and exit

Arguments:
  INPUT  Path to input file, multiple files separate by space [bam|sam|SJ.out.tab|gmap align|SJ]
```
唯有`-b`参数比较复杂，此处使用BAM/SAM文件来计算IR的PSI。因此，需要通过该参数指定BAM/SAM所在的路径。
- 若Input文件为BAM/SAM文件，无需指定该参数
- 若Input文件为STAR输出SJ.out.tab文件，如果没有做任何调整，则SJ.out.tab文件与BAM文件位于同一目录下，此时自动进行匹配，无需指定该参数。
- 直接指定BAM/SAM文件则直接使用该文件
- 指定BAM/SAM所在的文件夹路径，则根据input file的名称自动匹配对应的BAM/SAM文件可能存在一定误差(原本主要针对STAR的输出结果)


#### 3. Identify alternative splicing events based on SMRT-seq data
```bash
Usage: sms [OPTIONS] [INPUT]...

  Identify alternative splicing events from SMRT-seq

Options:
  -r, --reference PATH                  Path to reference file [gtf|gff3]
  -b, --bam PATH                         Path to BAM/SAM file
  -o, --output PATH                      Path to output file
  -p, --process INT                        Number of processes to use
  -e INT                                         The error to identify whether AS event exists
                               [default: 3bp]
  -c, --count INT                           Filter low abundance junctions [default: 0]
  --overlap-ref-reads FLOAT        Minimal overlap level to match reference with
                               reads [default: 90.0]
  --overlap-exon-intron FLOAT     Minimal overlap level between exon with intron
                              required for intron retention         	
                              identification [default: 90.0]
  -v, --verbose                                Enable detailed messages
  -h, --help                                     Show this message and exit

Arguments:
  INPUT  Path to input file, multiple files separate by space [BAM|SAM|gmap align|SJ]
```
基本参数与sgs模式雷同

`-b`若输入为BAM/SAM文件则无需指定该参数。

#### 4. Construct isoforms based on SMRT-seq data
```bash
Usage: iso [OPTIONS] [INPUT]...

  Construct Isoforms through SMRT-seq data

Options:
  -r, --reference PATH       	     Path to reference file [gtf|gff3]
  -o, --output PATH          	      Path to output file
  --overlap-ref-reads FLOAT  	Minimal overlap level to match reference with reads 
						  [default: 90.0]
  -v, --verbose             		 Enable detailed messages
  -h, --help                 		    Show this message and exit

Arguments:
  INPUT  Path to input BAM/SAM files, multiple files separate by space
```
## FAQ
1. If there are too many input files, may result in failed to check input file format, and so on.
2. subtypes of SE events
  - exact: the junciton just cross exact one exons -> junciont : 100-200, then the exon is 101-199
  - part: only one side of junciton match with known exon
  - other: junctions cross only one exon, but both start and end sites do not match with junctions
  - single: junction cross multiple exons of same transcript
  - multi: junction cross multiple exons from different transcripts
  
