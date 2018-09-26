package dsu.second.index

import dsu.carrier.Exons
import dsu.progressbar.ProgressBar
import java.io.File
import java.io.IOException
import java.util.*



/**
 * @author Zhangyiming
 * @since 2018.09.03
 * @version 20180903
 * 从Gtf文件中提取所有的exon
 */

class GffIndex(infile: String) : AnnotationIndex (infile) {

    override fun getSource(info: List<String>): Map<String, String> {
        val results: MutableMap<String, String> = mutableMapOf()

        for (inf in info[0].split(";")) {
            val tmp = inf.split("=")

            if ( ":" in tmp[1] ) {
                results[tmp[0]] = tmp[1].split(":")[1]
            } else {
                results[tmp[0]] = tmp[1]
            }
        }

        return results
    }

    override fun readExon() {

        try {
            val reader = Scanner(File(this.infile))

            val pb = ProgressBar(message = "Reading exon from Gff")
            val geneTranscript = mutableMapOf<String, String>()
            while (reader.hasNext()) {
                val line = reader.nextLine()

                if ( line.startsWith("#") ) {
                    continue
                }

                val lines = line.split("\\s+".toRegex())
                pb.step()

                val sources = this.getSource(lines.subList(8, lines.size))
                if ( "transcript_id" in sources.keys ) {
                    geneTranscript[sources["ID"]!!] = sources["Parent"]!!
                } else if (lines[2] == "exon") {
                    val tmpExon = Exons(
                            chromosome = lines[0],
                            start = lines[3].toInt(),
                            end = lines[4].toInt(),
                            strand = lines[6].toCharArray()[0],
                            exonId = sources["exon_id"]!!
                    )

                    tmpExon.source["transcript"] = sources["Parent"]!!
                    tmpExon.source["gene"] = geneTranscript[sources["Parent"]]!!

                    val tmp = mutableListOf(tmpExon)
                    val key = "${lines[0]}${lines[6]}"
                    if ( this.data.containsKey(key) ) {
                        this.data[key]!!.addAll(tmp)
                    } else {
                        this.data[key] = tmp
                    }
                }
            }

            reader.close()
            pb.close()
        }catch (e: IOException) {
            this.logger.error(e.toString())
        }

    }
}