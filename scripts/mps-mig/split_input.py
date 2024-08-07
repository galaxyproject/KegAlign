#!/usr/bin/env python3

import os
import argparse
import subprocess
import heapq

def chunk_file(
    file_name,
    output_dir,
    chunk_size=int((10 * 10**6) / 50),
    overlap_size=int((10 * 10**3) / 50),
    end_chr=None,
    save_split=False,
    single_chr=False
):
    file = open(file_name)
    data = file.readlines()  # each line (in .fa file) has up to 50 characters
    file.close()

    title = data[0]

    if single_chr:
        for line in reversed(data):
            if ">" in line:
                data.remove(line)
        data.insert(0, title)

    assert (">" in title)
    if end_chr is not None:
        split_index = -1
        for ind, line in enumerate(data):
            if end_chr == line.strip():
                split_index = ind
                break
        if split_index == -1:
            print(f"Stopping; chromosome {end_chr} not found")
            exit(0)
        print(f"splitting data to line {split_index} out of {len(data)}")
        data = data[:split_index]
        if save_split:
            save_split_file = file_name + f".to_{end_chr.strip()[1:]}"
            print(f"Saving split file {file_name} to {save_split_file}")
            with open(save_split_file, "w") as f:
                f.write(data)

    for line_index in range(0, len(data), chunk_size - overlap_size):
        end_line_index = min(len(data), line_index + chunk_size)
        block_file_name = os.path.join(
            output_dir, f"block_part{line_index}-{end_line_index}"
        )

        with open(block_file_name, "w") as out_file:
            to_write = data[line_index:end_line_index]

            if ">" not in to_write[0]:
                out_file.write(title)
            out_file.writelines(to_write)

            #  find last chromosome name for next file
            for line in reversed(to_write[:-overlap_size]):
                if ">" in line:
                    title = line
                    break

def convert_to_2bit(dir_name):
    process_list = []
    dir_files = os.listdir(dir_name)
    
    for file in dir_files:
        new_name = file + '.2bit'
        command = f"faToTwoBit {os.path.join(dir_name, file)} {os.path.join(dir_name, new_name)}"# && rm {os.path.join(dir_name, file)}"
        process = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True, shell = True)
        process_list.append(process)
        
    for process in process_list:
        stdout, stderr = process.communicate()
        print(stdout, end='')
        print(stderr, end='')
    
def split_chr(
    file_name,
    output_dir,
    num_chunks,
    write_to_output_dir = True,
    verbose = True
):
    
    file = open(file_name, 'r')
    data = file.read()  # each line (in .fa file) has up to 50 characters
    file.close()

    chrs = data.split('>')[1:]
    
    # Longest-processing-time-first Algorithm

    chrs = sorted(chrs, key=lambda x: (-len(x), x)) # sort ascending order
    
    chunk_size_list = []
    
    pq = []
    for i in range(num_chunks):
        heapq.heappush(pq, (0, i)) # bin size and bin id
    
    files = {}
    for i in range(num_chunks):
        files[i] = []
    
    for i, chr_ in enumerate(chrs):
        # get smallest file
        size, bin = heapq.heappop(pq)
        size += len(chr_)
        heapq.heappush(pq, (size, bin))
        files[bin].append(i)
        seen_inds = [] # sanity checking

    seen_inds = [] # for sanity checking
    for bin_id, chr_indexes in files.items():
        to_write = ''
        bin_size = 0
        for ind in chr_indexes:
            to_write += '>' + chrs[ind] + '\n'
            assert(ind not in seen_inds) # sanity check
            seen_inds.append(ind)
            bin_size += len(chrs[ind])
        if verbose:
            print(f"chunk_{bin_id} num bp {bin_size}")
        chunk_size_list.append(bin_size)
        block_file_name = os.path.join(output_dir, f'chunk_{bin_id}')
        if write_to_output_dir:
            with open(block_file_name, "w") as out_file:
                out_file.write(to_write)
    if verbose:
        print(f"packed {len(chrs)} chromosomes into {num_chunks} bins")
    assert(len(seen_inds) == len(chrs))
    assert(len(chunk_size_list) == num_chunks)
    return(chunk_size_list)

def parallel_wrapper(pass_list):
    return split_chr(*pass_list)

def mse(data, base):
    mse = [(base-i)**2 for i in data]
    return sum(mse)/len(mse)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str)
    parser.add_argument(
        "--out", type=str, default=None
    )
    parser.add_argument(
        "--end_chr", type=str, default=None, help="Use data up until this chromosome"
    )
    parser.add_argument(
        "--save_split_inputs", type=bool, default=False, help="Save split input files"
    )
    parser.add_argument(
        "--single_chr", type=bool, default=False, help="Use only first chromosome name for all chromosomes"
    )
    parser.add_argument(
        "--to_2bit", type=bool, default=False, help="Convert data into .2bit format"
    )
    parser.add_argument(
        "--max_chunks", type=str, default=None, help="Calculate chunk size based on input size"
    )
    parser.add_argument(
        "--split_chr", type=bool, default=True, help="Make each chromosome a separate partition"
    )
    parser.add_argument(
        "--goal_bp", type=int, default=False, help="goal basepairs across chunks, calculate best using MSE. Check up to --max_chunks number of bins"
    )
    args = parser.parse_args()

    target_file = args.input #'/scratch/mdl/WGA_tests/chr4.fa' #'/scratch/mdl/WGA_tests/galGal3.chr4.fa'

    target_block_dir = args.out  #'/scratch/mdl/WGA_tests/target_blocked_chr4/'
    
    query_chunk_size = (10 * 10**6) / 50
    target_chunk_size = (10 * 10**6) / 50
    overlap_size=int((10 * 10**3) / 50)

    if not os.path.exists(target_block_dir):
        print(f"making dir for {target_block_dir}")
        os.makedirs(target_block_dir)

    if args.max_chunks and not args.split_chr:
        import math
        chunk_size = args.max_chunks.split(',')
        if len(chunk_size) == 1:
            num_chunks_target = int(chunk_size[0])
        else:
            num_chunks_target = int(chunk_size[0])

        with open(target_file, 'r') as f:
            target_lines = len(f.readlines())
        
        def calculate_chunk_size(N, num_chunks, overlap):
            return math.ceil(N / num_chunks) + overlap
          
        target_chunk_size = calculate_chunk_size(target_lines, num_chunks_target, overlap_size) # each line has 50 base pairs
        print(f"overlap = {overlap_size}")
        print(f"target num lines = {target_lines}")
        print(f"target chunk size = {target_chunk_size}")

    if args.end_chr is None:
        chr_end_target = None
    else:
        chr_end_target = args.end_chr
    
    if args.split_chr:
        from math import inf
        bin_count = int(args.max_chunks)
        if args.goal_bp:
            USE_MULTIPROCESSING = True
            best_bin_count = None
            if USE_MULTIPROCESSING:
                import multiprocessing as mp
                pass_list = [(target_file, '', i, False, False) for i in range(1, bin_count+1)]
                bins_list = []
                with mp.Pool() as p:
                    bins_list = p.map(parallel_wrapper, pass_list)
                    
                prev_avg = 0
                goal_bp = int(args.goal_bp)
                best_bin_count = -1
                best_bin_loss = inf
                for i, bins in enumerate(bins_list):
                    i+=1
                    #avg = sum(bins)/len(bins)
                    #loss = abs(goal_bp - avg)
                    loss = mse(bins, goal_bp)
                    print(f"* bin count {i}, mse {int(loss)}, bins {bins}")
                    if loss < best_bin_loss:
                        best_bin_count = i
                        best_bin_loss = loss
                    
            else: 
                prev_avg = 0
                goal_bp = int(args.goal_bp)
                best_bin_count = -1
                best_bin_loss = inf
                for i in range(1, bin_count+1):
                    bins = split_chr(target_file, target_block_dir, num_chunks=i, write_to_output_dir=False, verbose=False)
                    avg = sum(bins)/len(bins)
                    loss = abs(goal_bp - avg)
                    print(f"bin count {i}, avg difference {int(loss)}")# (avg_bp {int(avg_bp)}, current_avg_bp {avg}, bins {bins})")
                    if loss < best_bin_loss:
                        best_bin_count = i
                        best_bin_loss = loss

            bin_count = best_bin_count
            
        print(f'bin_count = {bin_count}, loss={best_bin_loss}')
        split_chr(
            target_file,
            target_block_dir,
            bin_count
        )
    
    else:
        chunk_file(
            target_file,
            target_block_dir,
            end_chr=chr_end_target,
            save_split=args.save_split_inputs,
            single_chr=args.single_chr,
            chunk_size=target_chunk_size,
            overlap_size=overlap_size
        )

    if args.to_2bit:
        import subprocess
        # Convert to 2bit format
        convert_to_2bit(target_block_dir)
        
'''
python split_input.py --input /home/mdl/abg6029/WGA_tests/inputs/hg38.fa --out /home/mdl/abg6029/WGA_tests/inputs/blocked_20_chr_hg38/ --to_2bit true --split_chr true --max_chunks 20

python split_input.py --input hg38.fa --out ./blocked_20_chr_hg38 --to_2bit true --goal_bp 200000000 --max_chunks 20

'''