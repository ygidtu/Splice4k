package dsu.third.template


import org.apache.log4j.Logger
import java.io.File
import java.io.IOException
import java.io.PrintWriter
import dsu.carrier.SpliceEvent
import dsu.carrier.SpliceGraph
import dsu.progressbar.ProgressBar



/**
 * @since 2018.06.21
 * @version 20180918
 * @author Zhang Yiming
 *
 * 根据基因和Reads的配对情况，找出其中的可变剪接情况
 */


class SJFinder(
        private val template: GeneReadsCoupler
) {
    private val logger = Logger.getLogger(SJFinder::class.java.toString())
    private val results = hashMapOf<SpliceEvent, String>()

    init {
        this.identifySJ()
    }

    /**
     * 识别各种可变剪接类型
     */
    private fun identifySJ() {


        // 这个gap就是为了控制输出一个合适的进度条的
        val pb = ProgressBar(this.template.templates.size.toLong(), "Splice events identifying at")
        for (pair in this.template.templates) {

            pb.step()

            val graph = SpliceGraph(
                    chromosome = pair.template.chromosome,
                    strand = pair.template.strand
            )

            for (i in pair.getReadsExons() ) {
                for ( j in 0..(i.size - 1) step 2 ) {
                    graph.addEdge(start = i.sorted()[j], end = i.sorted()[j + 1])
                }

            }

            for ( i in graph.identifyAS(true) ) {
                this.results[i] = pair.template.transcriptId
            }
        }
        pb.close()
    }


    /**
     * 保存至文件
     * @param outfile 输出文件的路径
     * @return
     */
    fun saveTo(outfile: String) {
        val outFile = File(outfile).absoluteFile

        var writer = PrintWriter(System.out)
        try{
            if (!outFile.parentFile.exists()) outFile.parentFile.mkdirs()

            writer = PrintWriter(outFile)

            for ( (k, v) in this.results ) {
                writer.println("$k\t$v")
            }

        } catch (err: IOException) {
            logger.error(err.message)
            for (i in err.stackTrace) {
                logger.error(i)
            }
        } finally {
            writer.close()
        }
    }
}