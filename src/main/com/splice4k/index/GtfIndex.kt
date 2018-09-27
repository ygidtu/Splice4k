package splice4k.index

import splice4k.base.Exons
import splice4k.base.Genes
import splice4k.progressbar.ProgressBar
import java.io.File
import java.io.IOException
import java.util.*


/**
 * @author Zhangyiming
 * @since 2018.09.03
 * @version 20180903
 * 从Gtf文件中提取所有的exon
 */

class GtfIndex(
        infile: String,
        smrt: Boolean = true
) : AnnotationIndex ( infile, smrt ) {

    override fun getSource(info: List<String>): Map<String, String> {
        val results: MutableMap<String, String> = mutableMapOf()

        for (i in 0..(info.size - 1) step 2) {
            try{
                results[info[i]] = info[i+1]     // get rid of useless characters
                        .replace("\"", "")
                        .replace(";", "")
            } catch ( e: IndexOutOfBoundsException ) {
                continue
            }


        }
        return results
    }

    override fun readExon() {

        try {
            val reader = Scanner(File(this.infile))

            this.logger.info("Reading from ${this.infile}")
            val pb = ProgressBar(message = "Reading exon from Gtf")

            val exons: MutableMap<String, MutableList<Int>> = mutableMapOf()

            while (reader.hasNext()) {
                val line = reader.nextLine()
                val lines = line.split("\\s+".toRegex())
                pb.step()

                if ( line.startsWith("#") || lines[0] !in arrayOf("exon", "transcript") ) {
                    continue
                }

                val sources = this.getSource(lines.subList(8, lines.size))

                if ( this.smrt && lines[2] == "transcript" ) {
                    this.transcripts.add(Genes(
                            chromosome = lines[0],
                            start = lines[3].toInt(),
                            end = lines[4].toInt(),
                            strand = lines[6].toCharArray()[0],
                            information = sources
                    ))
                }

                if (lines[2] == "exon") {

                    val tmp = Exons(
                            chromosome = lines[0],
                            start = lines[3].toInt(),
                            end = lines[4].toInt(),
                            strand = lines[6].toCharArray()[0],
                            exonId = sources["exon_id"]!!
                    )

                    tmp.source["transcript"] = sources["transcript_id"]!!
                    tmp.source["gene"] = sources["gene_id"]!!

                    val tmpExons = mutableListOf(tmp)
                    val key = "${lines[0]}${lines[6]}"
                    if ( this.data.containsKey(key) ) {
                        this.data[key]!!.addAll(tmpExons)
                    } else {
                        this.data[key] = tmpExons
                    }

                    if ( this.smrt ) {
                        val tmpExon = mutableListOf(tmp.start, tmp.end)

                        if (exons.containsKey(sources["Parent"]!!)) {
                            tmpExon.addAll(exons[sources["Parent"]!!]!!)
                        }
                        exons[sources["Parent"]!!] = tmpExon
                    }
                }
            }

            reader.close()
            pb.close()


            if ( this.smrt ) {
                for (i in transcripts) {
                    if (exons.containsKey(i.transcriptId)) {
                        i.exons.addAll(exons[i.transcriptId]!!)
                    }
                }
            }

        }catch (e: IOException) {
            this.logger.error(e.toString())
        }

    }
}