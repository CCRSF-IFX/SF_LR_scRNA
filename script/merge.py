import pandas as pd
import argparse
import logging

def read_classification_file(input_file):
    logging.info(f'Reading classification file: {input_file}')
    return pd.read_csv(input_file, sep='\t')

def read_count_file(count_file):
    logging.info(f'Reading count file: {count_file}')
    return pd.read_csv(count_file, sep='\t')

def read_isoform_map_file(map_file):
    logging.info(f'Reading isoform map file: {map_file}')
    return pd.read_csv(map_file, sep='\t', header=None, names=['isoform', 'read_ids'])

def read_read_tags_file(read_tags_file):
    logging.info(f'Reading read tags file: {read_tags_file}')
    return pd.read_csv(read_tags_file, sep='\t')

def replace_ids(read_ids, read_id_to_barcode):
    ids = read_ids.split(',')
    corrected_barcodes = [read_id_to_barcode.get(read_id.strip()) for read_id in ids if read_id.strip() in read_id_to_barcode]
    return ','.join(corrected_barcodes)

def replace_read_ids_with_barcodes(isoform_map_df, read_tags_df):
    logging.info('Replacing read IDs with corrected barcodes')
    read_tags_df['trimmed_read_id'] = read_tags_df['read_id'].str.replace(r'_\d+$', '', regex=True)
    read_id_to_barcode = dict(zip(read_tags_df['trimmed_read_id'], read_tags_df['corrected_barcode']))
    isoform_map_df['corrected_read_ids'] = isoform_map_df['read_ids'].apply(replace_ids, read_id_to_barcode=read_id_to_barcode)
    isoform_map_df = isoform_map_df.drop(columns=['read_ids'])
    return isoform_map_df

def merge_dataframes(classification_df, count_df, isoform_map_df):
    logging.info('Merging count and isoform map dataframes')
    merged_count_isoform_df = pd.merge(count_df, isoform_map_df, left_on='ids', right_on='isoform', how='left')
    merged_count_isoform_df = merged_count_isoform_df.drop(columns=['isoform'])
    merged_count_isoform_df.columns.values[-2] = 'FL counts'
    merged_count_isoform_df.columns.values[-1] = 'barcodes'
    merged_count_isoform_df['isoform_extract'] = merged_count_isoform_df['ids'].str.split('_').str[:2].str.join('_')
    
    logging.info('Merging classification dataframe with merged count and isoform map dataframe')
    merged_df = pd.merge(classification_df, merged_count_isoform_df, left_on='isoform', right_on='isoform_extract', how='left')
    merged_df = merged_df.drop(columns=['ids', 'isoform_extract'])
    return merged_df

def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    parser = argparse.ArgumentParser(description='Merge SQANTI classification, count, and isoform map files.')
    parser.add_argument('-i', '--input', required=True, help='Path to the SQANTI classification file')
    parser.add_argument('-c', '--count', required=True, help='Path to the count file')
    parser.add_argument('-m', '--map', required=True, help='Path to the isoform.read.map.txt file')
    parser.add_argument('-r', '--read-tags', required=True, help='Path to the read_tags.tsv file')
    parser.add_argument('-o', '--output', required=True, help='Path to the output file')
    
    args = parser.parse_args()

    logging.info('Starting the merging process')
    classification_df = read_classification_file(args.input)
    count_df = read_count_file(args.count)
    isoform_map_df = read_isoform_map_file(args.map)
    read_tags_df = read_read_tags_file(args.read_tags)

    logging.info('Replacing read IDs with corrected barcodes')
    isoform_map_df = replace_read_ids_with_barcodes(isoform_map_df, read_tags_df)
    
    logging.info('Merging all dataframes')
    merged_df = merge_dataframes(classification_df, count_df, isoform_map_df)
    
    logging.info(f'Writing output to file: {args.output}')
    merged_df.to_csv(args.output, sep='\t', index=False)
    logging.info('Merging process completed successfully')

if __name__ == '__main__':
    main()