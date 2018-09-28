package com.splice4k.index


import com.splice4k.base.Genes
import com.splice4k.base.SpliceGraph
import com.splice4k.progressbar.ProgressBar
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import java.io.File
import java.io.IOException
import java.io.PrintWriter


/**
  * @author Zhangyiming
  * @since 2018.09.03
  * @version 20180903
  * Extract Splice Junctions from Bam/Sam file
  */


class BamIndex(
        infile: String,
        val silent: Boolean = true,
        val unique: Int = 1,
        val smrt: Boolean = false,
        filter: Int
): SJIndex(infile = infile, filter = filter, star = false) {
    val transcripts = mutableListOf<Genes>()

    /**
     * private function
     * 从单条bam中提取cigar中的N区域
     * @param record 单条SAM/BAM的信息
     * @return 列表，记录了所有的intron的边界信息
     */
    private fun extractSpliceFromCigar( record: SAMRecord ): List<Int> {
        val results: MutableList<Int> = mutableListOf()
        var position = record.alignmentStart
        val tmp = mutableListOf<Char>()

        for (i in record.cigar.toString()) {
            if (i in '0'..'9') {  // 如果是数字，就加到list中
                tmp.add(i)
            } else {
                if (tmp.size == 0) {
                    continue
                }

                // Soft clip以及insertion的两种区域都不记载在alignment之内
                if (i != 'S' && i != 'I') {
                    position += tmp.joinToString(prefix = "", postfix = "", separator = "").toInt()
                }

                if (i == 'N') {
                    results.add(
                            position - tmp.joinToString(prefix = "", postfix = "", separator = "").toInt()
                    )

                    results.add(position - 1)
                }
                tmp.clear()
            }
        }

        return results
    }


    /**
     * 获取bam文件中的所有SJ及其freq数值
     * 在处理SMRT数据时，顺便构建transcript及外显子列表
     * @return Map<GenomicLoci, Int>
     */
    override fun getAllSJ() {
        val tmpReader =  SamReaderFactory
                .makeDefault()
                .open(File(this.infile))
                .iterator()

        this.logger.info("Reading from ${this.infile}")
        val pb = ProgressBar(message = "Reading from Bam")
        for ( record in tmpReader) {

            pb.step()

            val spliceSites = this.extractSpliceFromCigar(record)

            if (spliceSites.isEmpty()) {
                continue
            }

            val strand = when(record.readNegativeStrandFlag) {
                true -> '-'
                false -> '+'
            }

            // SGS构建junctions map
            val tmpGraph = when( this.data.containsKey("${record.referenceName}$strand")  ) {
                true -> this.data["${record.referenceName}$strand"]!!
                else -> SpliceGraph(record.referenceName, strand)
            }

            for ( i in 0..(spliceSites.size - 1) step 2) {
                tmpGraph.addEdge(start = spliceSites[i], end = spliceSites[i + 1])
            }

            this.data["${record.referenceName}$strand"] = tmpGraph

            // SMRT 构建list of transcripts
            if ( smrt ) {
                // 判断reads是否为unique mapped
                if (record.hasAttribute("NH")) {

                    val mapped = record.getAttribute("NH").toString().toInt()

                    if ( mapped <= 0 || mapped > this.unique) continue
                } else {
                    // 没有NH标签的reads，通常也会造成其他错误，因此直接放弃
                    if (!this.silent) this.logger.warn("${record.readName} does not have attribute NH")
                    continue
                }


                val junctions = this.extractSpliceFromCigar(record)
                // init Genes
                val tmpGene = Genes(
                        chromosome = record.referenceName,
                        start = record.alignmentStart,
                        end = record.alignmentEnd,
                        geneName = record.readName,
                        strand = strand
                )
                // construct exon to transcripts
                tmpGene.exons.addAll(junctions)

                this.transcripts.add(tmpGene)
            }
        }
        pb.close()

        // filtering
        if ( this.filter > 0 ) {
            this.logger.info("Filtering splice junctions")
            this.data.values.forEach{ it.filter(this.filter) }
        }
    }


    /**
     * 将提取到的Junctions输出到文件
     * @param output 输入文件
     */
    fun writeTo( output: File ) {

        try{
            if (!output.absoluteFile.parentFile.exists()) {
                output.absoluteFile.parentFile.mkdirs()
            }

            val writer = PrintWriter(output)

            for ( v in this.data.values ) {
                writer.print(v.toString())
            }
            writer.close()
        } catch (err: IOException) {
            logger.error(err.message)
            for (i in err.stackTrace) {
                logger.error(i)
            }
        }
    }

 }