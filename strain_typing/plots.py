import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pylab
from matplotlib.patches import Rectangle
matplotlib.rcParams['lines.linewidth'] = 1


def produce_histograms(histo, coverage, name, cutoff, cutoff_count, barcode, estimate_genome_size, output_path):
    """

    :param histo:
    :param coverage:
    :param name:
    :param cutoff:
    :param cutoff_count:
    :param barcode:
    :param estimate_genome_size:
    :param output_path:
    :return:
    """
    fig = pylab.figure()
    fig.set_size_inches(10, 3)

    ax = fig.add_subplot(1, 1, 1)

    freq, count = [], []
    histo = dict(histo)
    for i in range(1, max(histo)):
        freq.append(i)
        if i in histo:
            count.append(histo[i])
        else:
            count.append(0)
    ylim = None
    try:
        ylim = int(histo[coverage] * 2.5) + 1
    except KeyError:
        ylim = int(max(count[int(cutoff_count):]) * 2.5) + 1

    xlim = coverage * 2

    width = 0.8  # the width of the bars

    ax.set_title("{0}\n{1}".format(name, barcode), fontsize=13)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    ax.bar(freq, count, width, color='g', linewidth=.5,)
    ax.set_xlim(0, xlim)
    ax.set_ylim(0, ylim)
    # cutoff
    ax.plot([cutoff, cutoff], [0, ylim], "--", color='#C24641')  # add threshold line
    ax.annotate('<{0} excluded (count at exclusion {1:,}'.format(cutoff, cutoff_count),
                xy=(cutoff, ylim * .95),
                xytext=(cutoff + .4, ylim * .95),
                fontsize=6, color='#C24641', rotation=90)

    ax.add_patch(Rectangle((0, 0), cutoff, ylim, alpha=.5, facecolor='#C24641', linewidth=0))

    # estimated coverage
    ax.plot([coverage, coverage], [0, ylim], "k--")  # add threshold line
    ax.annotate('{:.1f}x estimated coverage'.format(coverage), xy=(coverage, ylim * .90),
                xytext=(coverage + .8, ylim * .90), fontsize=8, color='#483C32')

    # estimate genome size
    ax.annotate('Estimated genome size\n{:,} bp'.format(estimate_genome_size),
                xy=(xlim * .90, ylim * .80),
                xytext=(xlim * .90, ylim * .60), fontsize=8, color='#483C32', horizontalalignment='right')

    pylab.xticks(fontsize=9)
    pylab.yticks(fontsize=9)

    ax.set_xlabel("kmer frequency", fontsize=12)
    ax.set_ylabel("kmer count", fontsize=12)
    try:
        fig.savefig(output_path, format='png', bbox_inches="tight", dpi=400)
    except ValueError:
        try:
            fig.savefig(output_path, format='png', bbox_inches="tight", dpi=100)
        except ValueError:
            output_path = output_path.replace(".png", ".pdf")
            fig.savefig(output_path, format='pdf', bbox_inches="tight")
            output_path = None
    pylab.close()
    return output_path

