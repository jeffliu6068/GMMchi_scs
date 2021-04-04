from matplotlib import pyplot as plt
import numpy as np
import scprep
from tqdm import tqdm_notebook as tqdm
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats, optimize, interpolate
import scipy


def doublet_empty_filter(input_data, doublet_cutoff=25000, empty_cutoff=500):
    print('removing barcodes that are highly likely to be doublets or empty...')

    # remove counts/UMI too many or too little
    f, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))

    original_len = input_data.shape[0]
    input_data = scprep.filter.filter_library_size(
        input_data, cutoff=empty_cutoff, keep_cells='above')
    print('{:.2%} barcodes were likely unused or empty...'.format(
        (original_len-len(input_data))/original_len))

    # before filtering
    sum_counts = input_data.sum(axis=1)
    ax1.hist(np.log2(sum_counts+1), 100, alpha=0.5, histtype='bar', ec='black')
    ax1.set_xlabel('Library size')
    ax1.set_ylabel('# of cells')
    ax1.set_title('Before filtering')

    # remove those with too many counts/UMI because may be dublets or too little <500
    input_data = scprep.filter.filter_library_size(
        input_data, cutoff=doublet_cutoff, keep_cells='below')

    print('{} barcodes left...'.format(len(input_data)))

    # after filtering
    sum_counts = input_data.sum(axis=1)
    ax2.hist(np.log2(sum_counts+1), 100, alpha=0.5, histtype='bar', ec='black')
    ax2.set_xlabel('Library size')
    ax2.set_ylabel('# of cells')
    ax2.set_title('After filtering')

    plt.tight_layout()
    plt.show()

    return input_data


def library_complexity_filter(input_data):

    counts = 0
    print('measuring library complexity for each barcode...')
    complexity_full = input_data[input_data > counts].count(axis=1)
    complexity_full = np.sort(complexity_full.values)

    complexity_full_lg = np.log2(complexity_full)
    plt.hist(complexity_full_lg, 40, alpha=0.5, histtype='bar', ec='black')
    plt.title(
        'Distribution of the number of UMIs a barcode has with counts >{}'.format(counts))
    plt.show()

    # initiate
    complexity_full_lg_ = complexity_full_lg.copy()

    skewness_full = []

    print('removing barcodes with low library complexity...')
    for i in tqdm(complexity_full_lg_):
        complexity_full_lg_ = complexity_full_lg_[1:]

        skewness_full.append(scipy.stats.skew(
            complexity_full_lg_, bias=True, nan_policy='propagate'))

    plt.plot(np.arange(0, len(skewness_full)), skewness_full)
    plt.title('Skewness of distribution')
    plt.show()

    # check if the skew is negative, if not skip
    if skewness_full[0] < 0:
        product_full = [skewness_full[i]*skewness_full[i+1]
                        for i in range(0, len(skewness_full)-1)]
        intercept = [index for index, i in enumerate(product_full) if i < 0]

        print('Intercept of skewness crossing 0 is: {} or 2^{}'.format(
            intercept[0], round(complexity_full_lg[intercept[0]], 2)))

        threshold = complexity_full_lg[intercept[0]]

        complexity_full = input_data[input_data > counts].count(axis=1)

        a = pd.DataFrame(complexity_full, index=input_data.index,
                         columns=['complexity'])
        input_data = input_data.loc[a[a['complexity'] > 2**threshold].index]

    print('{} barcodes left...'.format(len(input_data)))

    return input_data


def housekeepinggene_filter(input_data):
    # housekeeipng genes that are most abundant
    from tqdm import tqdm_notebook as tqdm
    from matplotlib import pyplot as plt
    import numpy as np
    import GMMchi

    print('removing barcodes with housekeeping genes under cutoff...')

    # get first 10 most abundant housekeeping genes that is not mitochondrial gene

    housekeeping = ['GAPDH', 'ACTG1', 'ACTB']
    housekeeping_df = input_data[housekeeping].T

    housekeeping_cutoff = {}

    for gene in tqdm(housekeeping):
        info, classif, categories, chi, bins, f = GMMchi.GMM_modelingt(gene, housekeeping_df, log2transform=True,
                                                                filt=0, meanf=0, stdf=0,
                                                                graphs=True, verbosity=False, chisquaremethod=True, Single_tail_validation=False)
        try:
            housekeeping_cutoff[gene] = info[2][0]
        except:
            housekeeping_cutoff[gene] = info[2]

    # transform the dictionary of cutoffs to 2** and also remove those that are unimodal
    housekeeping_cutoff = {k: 2**v for k,
                           v in housekeeping_cutoff.items() if v is not None}

    # transpose housekeeping_df
    housekeeping_df = housekeeping_df.T

    # remove housekeeping genes that is below cutoff
    for k, v in housekeeping_cutoff.items():
        housekeeping_df = housekeeping_df[housekeeping_df[k] > v]

    # get the index of the single cells UMI that is left after filtering
    filtered_housekeeping_index = housekeeping_df.index

    input_data = input_data.loc[filtered_housekeeping_index]

    print('{} barcodes left...'.format(len(input_data)))

    return input_data


def mitochondrial_gene_filter(input_data, mito_percentile=95):
    print('removing barcodes with mitochondrial genes that are in the top 5 percentile...')

    # get mitochondrial genes
    mitochondrial_gene_list = np.array(
        [g.startswith('MT-') for g in input_data.columns])

    # get expression
    mito_exp = input_data.loc[:, mitochondrial_gene_list].mean(axis=1)

    # plotting
    fig, ax = plt.subplots(1, figsize=(6, 5))

    ax.hist(mito_exp, bins=100, alpha=0.5, histtype='bar', ec='black')
    ax.axvline(np.percentile(mito_exp, 95))
    ax.set_xlabel('Mean mitochondrial expression')
    ax.set_ylabel('# of cells')
    ax.set_title('Mitochondrial expression')

    fig.tight_layout()
    plt.show()

    # filter cells above 95th percentile
    input_data = scprep.filter.filter_values(
        input_data, values=mito_exp, percentile=mito_percentile, keep_cells='below')

    print('{} barcodes left...'.format(len(input_data)))

    return input_data


def housekeeping_normalization(input_data):
    housekeeping_sum = input_data[[
        'GAPDH', 'ACTB', 'ACTG1']].mean(axis=1).values

    return input_data.div(housekeeping_sum, axis=0)


def GMMchi_scs_pipeline(input_data, doublet_cutoff=25000, empty_cutoff=500, mito_percentile=95, normalization=True):
    print('GMMchi_scs_pipeline initiating...')
    # remove rare genes, if a gene is not expressed in at least 5 cells remove
    input_data = scprep.filter.remove_rare_genes(
        input_data, cutoff=0, min_cells=5)

    input_data = doublet_empty_filter(input_data, doublet_cutoff, empty_cutoff)
    input_data = library_complexity_filter(input_data)
    input_data = housekeepinggene_filter(input_data)
    input_data = mitochondrial_gene_filter(input_data, mito_percentile)

    if normalization == True:
        print('normalizing input data...')
        # normalization, # of UMIs / cell
        input_data = scprep.normalize.library_size_normalize(input_data)
        # normalization usingonly sum of housekeeping gene
        # input_data = housekeeping_normalization(input_data)

    print('GMMchi_scs_pipeline finished. Done!')
    return input_data


def UMAP_graph(input_data, verbose=True):
    # Dimension reduction and clustering libraries
    import umap.umap_ as umap
    import hdbscan
    import sklearn.cluster as cluster
    from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
    import time

    # UMAP

    time_start = time.time()
    standard_embedding = umap.UMAP(
        random_state=42, verbose=verbose).fit_transform(input_data)

    print('UMAP done! Time elapsed: {} seconds'.format(time.time()-time_start))

    df_subset = pd.DataFrame()

    df_subset['UMAP-2d-one'] = standard_embedding[:, 0]
    df_subset['UMAP-2d-two'] = standard_embedding[:, 1]

    df_subset.index = input_data.index

    return df_subset


def Label_graph(input_data, df_subset, label_list=[], thresholds=[], use_gene=True, boolean_visualization=True):
    # generate different hue for gene_list

    import seaborn as sns
    from matplotlib import pyplot as plt
    import numpy as np
    from matplotlib import colors

    plt.figure(figsize=(16, 10))

    if boolean_visualization == True:
        if use_gene == True:
            gene_interest = label_list

            hue = ['None']*len(input_data)
            # if threshold_list_boolean == True:
            # preset threshold as 0
            if thresholds == []:
                thresholds = [0]*len(label_list)

            for y, threshold_ in zip(gene_interest, thresholds):
                for index, x in enumerate(input_data[y]):
                    if x > threshold_:
                        if hue[index] != 'None':
                            hue[index] = hue[index] + '/{}+'.format(y)
                        else:
                            hue[index] = '{}+'.format(y)

            # else:
            #     for y in gene_interest:
            #         for index, (x, threshold) in enumerate(input_data[y], thresholds):
            #             if x>threshold:
            #                 if hue[index] != 'None':
            #                     hue[index] = hue[index]+ '/{}+'.format(y)
            #                 else:
            #                     hue[index] = '{}+'.format(y)

            key = pd.value_counts(hue).index

            color_dict = {x: y for x, y in zip(
                key, [x for x in sns.color_palette("tab10", len(key))])}
            color_dict['None'] = 'lightgrey'
            # color_dict['{}+'.format(label_list[0])] = 'orange'
            # color_dict['{}+'.format(label_list[1])] = 'steelblue'
            # color_dict['{}+/{}+'.format(label_list[0], label_list[1])] = 'red'

            # plot it out
            df_subset['hue'] = hue
        else:
            label_list = list(label_list)
            df_subset['hue'] = label_list
            key = pd.value_counts(label_list).index

            color_dict = {x: y for x, y in zip(
                key, [x for x in sns.color_palette("tab10", len(key))])}
            color_dict[0] = 'lightgrey'

        cmap = color_dict
        legend = 'full'

        aa = sns.scatterplot(
            x="UMAP-2d-one", y="UMAP-2d-two",
            palette=cmap,
            hue='hue',
            data=df_subset,
            legend=legend,
            alpha=1
        )
        plt.show()
    else:
        for gene in label_list:
            plt.figure(figsize=(16, 10))

            label_name = '{} log2(expression level)'.format(gene)

            df_subset[label_name] = np.round(np.log2(input_data[gene]+1), 2)

            # hence the first value 0, second value 0.01
            c = np.linspace(0, 1, 101)
            # For those values we store the colors from the "YlOrRd" map in an array
            color = plt.get_cmap('coolwarm', 101)(c)
            # We replace the first row of that array, by white
            color[0, :] = colors.to_rgba('lightgrey')
            # We create a new colormap with the colors
            cmap = colors.ListedColormap(color)

            legend = 'brief'

            sns.scatterplot(
                x="UMAP-2d-one", y="UMAP-2d-two",
                palette=cmap,
                hue=label_name,
                data=df_subset,
                legend=legend,
                alpha=1,
                hue_norm=(0, max(df_subset[label_name]))
            )
            plt.show()
