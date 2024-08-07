import subprocess
import os
import re
import argparse
import datetime
import time
from pynvml import *


class NamedPopen(subprocess.Popen):
     """
     Like subprocess.Popen, but returns an object with a .name member
     """
     def __init__(self, *args, name=None, **kwargs):
         self.name = name
         super().__init__(*args, **kwargs)

class GPU_queue():
    """
    """
    def __init__(self, device_names:list, uid_folder:str, uid_prefix='UID_', max_processes:list=None):
        self.queue = {}
        for i in device_names:
            self.queue[i] = []
        self.uid_folder = uid_folder
        self.submitted_uid = set()
        self.completed_uid_history = set()
        self.uid_prefix = uid_prefix
        self.max_processes = {}
        if max_processes != None:
            for i in range(len(device_names)):
                self.max_processes[device_names[i]] = max_processes[i]

    def submit(self, uid:str, device_name:str):
        self.queue[device_name].append(uid)
        self.submitted_uid.add(uid)

        if max_processes != None:
            if (self.max_processes[device_name] < len(self.queue[device_name])):
                print(f"WARNING for device {device_name}, max process = {self.max_processes[device_name]}, running process = {len(self.queue[device_name])}")
                assert(False)

    '''
    Removes processes that completed their GPU part but not necessarily CPU (LASTZ) part.
    '''
    def check_completion(self):
        completed_uid = set([f for f in os.listdir(self.uid_folder) if self.uid_prefix in f and os.path.isfile(os.path.join(self.uid_folder, f))])

        # check successfully completed jobs using file output from modified run_segalign script
        uids_in_progress = completed_uid - self.completed_uid_history
        self.remove_uids(uids_in_progress)

    def remove_uids(self, uid_list):
        for uid in uid_list:
            for device_queue in self.queue.values():
                if uid in reversed(device_queue):
                    device_queue.remove(uid)
                    self.completed_uid_history.add(uid)
                    continue

    def get_queue(self):
        return self.queue

    def get_running_uids(self):
        uids = []
        for device_queue in self.queue.values():
            for uid in device_queue:
                uids.append(uid)
        return uids

    def get_free_device_list(self):
        free_device = []
        for key, value in self.queue.items():
            running_process = len(value)
            max_processes = self.max_processes[key]
            if running_process < max_processes:
                free_device.append((key, running_process))
        return free_device


'''
Hold list of processes.
Removes complete processes to prevent too many file open errors.
'''
class Process_List():
    def __init__(self, max_processes = 256):
        self.processes = {} #key process.name (UID), value process
        self.stdout = ""
        self.stderr = ""
        self.max_processes = max_processes

    def append(self, process):
        print(f"== ADDING process {process.name}")
        prev_len = len(self.processes)
        assert(self.processes.get(process.name) is None)
        self.processes[process.name] = process


        if len(self.processes) >= self.max_processes:
            self.get_fails_and_check_completion()

    def _detect_failure(self, stderr):
        failures = ["core dumped", "Can't open", "cuda", "bad_alloc", "Aborted"]
        for f in failures:
            if f in stderr:
                if "cudaErrorCudartUnloading" in stderr:
                    return 'mem_err'
                else:
                    return 'other_err'
        return 'no_err'

    '''
    Check if processes are complete and saves outputs. Returns any process names (UID) that have failed
    '''
    def get_fails_and_check_completion(self):
        process_names = list(self.processes.keys())
        fails = []
        mem_fails = []
        for p_name in process_names:
            p = self.processes[p_name]
            if p.poll() is not None:
                print("\33[2K", end='') # clear line
                print(f"-- REMOVING process {p.name}")
                stdout, stderr = p.communicate() # get outputs. Blocks until process returns
                self.stdout += stdout
                self.stderr += stderr
                error_type = self._detect_failure(stderr)
                if error_type == 'other_err':
                    print(f"FAIL {p.name}")
                    fails.append(p.name)
                if error_type == 'mem_err':
                    print(f"MEM_FAIL {p.name}")
                    mem_fails.append(p.name)
                del self.processes[p.name]
        return fails, mem_fails

    '''
    Checks uids in uid_list and returns those that are complete
    '''
    def check_uid_completion(self, uid_list):
        complete_uids = []
        for uid_name in uid_list:
            process = self.processes.get(uid_name)
            assert(process.name == uid_name)
            if process is not None and process.poll() is not None:
                complete_uids.append(uid_name)
        return complete_uids
    
    '''
    Check if finished processes, by uid name, failed due to memory errors
    '''
    def check_uid_memfail(self, uid_list):
        mem_fails = []
        for uid_name in uid_list:
            process = self.processes.get(uid_name)
            assert(process.name == uid_name)
            assert(process is not None and process.poll() is not None)
            stdout, stderr = process.communicate()
            error_type = self._detect_failure(stderr)
            if error_type == 'mem_err':
                mem_fails.append(process.name)
        return mem_fails

    def get_output(self, terminate=False):
        output_stdout = ''
        output_stderr = ''
        for name, process in self.processes.items():
            if terminate:
                if process.poll() is None: # if process still running
                    process.terminate()
                    print(f"Error 1: Could not get results for process {process.name}")
                    stdout = ''
                    stderr = ''
                else:
                    try:
                        stdout, stderr = process.communicate()
                    except:
                        print(f"Error 2: Could not get results for process {process.name}")
                        stdout = ''
                        stderr = ''
            else:
                print(f"waiting for process {name}")
                try:
                    stdout, stderr = process.communicate() # get outputs. Blocks until process returns
                except:
                    print(f"Error 3: Could not get results for process {process.name}")
                    stdout = ''
                    stderr = ''
            self.stdout += stdout
            self.stderr += stderr
        return self.stdout, self.stderr

    def print_fails(self):
        print(f"FAILED: {self.fails()}")


def run_command(command):
    line_as_bytes = subprocess.check_output(f"{command}", shell=True)
    line = line_as_bytes.decode("ascii")
    return line

def get_nvidia_smi():
    return run_command('nvidia-smi')

def get_nvidia_smi_L():
    return run_command('nvidia-smi -L')

def get_processes(): # DEPRECATED nvidia-smi very slow to exec.
    x = get_nvidia_smi().split('Processes:')[1].split('|=======================================================================================|')[1].strip()
    processes = x.split('\n')[:-1]

    if 'No running processes found' in processes[0]:
        return []

    process_list = []
    for line in processes:
        s = line.split()
        gpu = s[1]
        gi_id = s[2]
        ci_id = s[3]
        name = s[6]
        process_list.append([gpu, gi_id, ci_id, name])
    return process_list


def get_processes_v2(use_MPS = False):
    deviceCount = nvmlDeviceGetCount()
    process_list = []
    for i in range(deviceCount):
        handle = nvmlDeviceGetHandleByIndex(i)
        if use_MPS:
            processes = nvmlDeviceGetMPSComputeRunningProcesses_v3(handle)
        else:
            processes = nvmlDeviceGetComputeRunningProcesses_v3(handle)
        for process in processes:
            process_list.append((i, process))
    return process_list



def get_uuids(): # TODO make this use nvml functions
    mig_uuids = re.findall(r'\(UUID: MIG(.*?)\)', get_nvidia_smi_L())
    mig_uuid_list = []
    for id in mig_uuids:
        mig_uuid_list.append(f"MIG"+id)
    return mig_uuid_list


def init_mig_dict(used_uuids): # TODO make this use nvml functions
    x = get_nvidia_smi().split('MIG devices:')[1]
    border_text = x.split('\n')[1].strip()

    x = x.split('=|')[1].strip()
    x = x.split('Processes:')[0].strip()
    #devices = x.split('+------------------+----------------------+-----------+-----------------------+')[:-1]
    devices = x.split(border_text)[:-1]
    mig_dict = {} # key = (gpu, gi_id, ci_id), value = mig_uuid
    mig_uuids = get_uuids()

    assert(len(devices) == len(mig_uuids))
    for i, device in enumerate(devices):
        s = device.split()
        gpu = s[1]
        assert(gpu.isdigit())
        gi_id = s[2]
        assert(gi_id.isdigit())
        ci_id = s[3]
        assert(ci_id.isdigit()) 
        mig_dev = s[4]
        assert(mig_dev.isdigit())
        uuid = mig_uuids[i]
        if uuid in used_uuids:
            mig_dict[(gpu, gi_id, ci_id)] = uuid
    return mig_dict

def get_free_mig(mig_dict): # return mig uuids
    process_list = get_processes()
    free_mig_list = list(mig_dict.values())

    for process in process_list:
        gpu, gi_id, ci_id, name = process
        mig_uuid = mig_dict.get((gpu, gi_id, ci_id))
        if mig_uuid != None:
            free_mig_list.remove(mig_uuid)
    return free_mig_list

def get_free_mig_v2(mig_dict, no_MIG=False):
    #print('==========')
    processes = get_processes_v2()
    free_mig_list = list(mig_dict.values())
    #for i, p in processes:
    #    print(f"gpu:{i}  {p}")
    for gpu_id, process in processes:
        if no_MIG:
            gi_id = 'N/A'
            ci_id = 'N/A'
        else:
            gi_id = process.gpuInstanceId
            ci_id = process.computeInstanceId
        mig_uuid = mig_dict.get((str(gpu_id), str(gi_id), str(ci_id)))
        if mig_uuid != None:
            #print(f"removing uuid {mig_uuid} with gpu_id:{gpu_id}, gi_id:{gi_id}, ci_id:{ci_id}  {free_mig_list}")
            try:
                free_mig_list.remove(mig_uuid)
            except:
                pass
                # TODO why does this happen?
    return free_mig_list

def get_mig_process_count(mig_dict, no_MIG=False, use_MPS=False): # similar to get_free_mig, but returns a dict with process count for each mig device
    #print('==========')
    processes = get_processes_v2(use_MPS)
    #if use_MPS:
    #    mig_process_count_dict = dict.fromkeys(mig_dict.values(), -1) # init -1 to ignore nvidia-cuda-mps process
    #else:
    mig_process_count_dict = dict.fromkeys(mig_dict.values(), 0)

    for gpu_id, process in processes:
        if no_MIG:
            gi_id = 'N/A'
            ci_id = 'N/A'
        else:
            gi_id = process.gpuInstanceId
            ci_id = process.computeInstanceId
        mig_uuid = mig_dict.get((str(gpu_id), str(gi_id), str(ci_id)))
        if mig_uuid != None:
            mig_process_count_dict[mig_uuid] += 1
    return mig_process_count_dict

def combine_results(output_dir, output_file, part_pattern='part_*.maf', remove=True):
    try:
        run_command(f"cat {os.path.join(output_dir,part_pattern)} > {output_file}")
        # "cat *.txt > all.txt"
        if remove:
            run_command(f"rm {os.path.join(output_dir,part_pattern)}")
    except:
        print("Could not combine results")

def remove_uids(uid_dir, uid_prefix):
    try:
        run_command(f"rm {os.path.join(uid_dir,uid_prefix)}*")
    except:
        print("Could not remove uids")

def init_mps(mig_list, mps_pipe_dir):
    for mig in mig_list:
        print(f"Initializing MPS server for {mig}")
        command = f'CUDA_VISIBLE_DEVICES={mig} CUDA_MPS_PIPE_DIRECTORY={os.path.join(mps_pipe_dir, mig)} nvidia-cuda-mps-control -d'
        run_command(command)

def destroy_mps(mig_list, mps_pipe_dir):
    for mig in mig_list:
        print(f"Destroying MPS server for {mig}")
        command = f'echo quit | CUDA_VISIBLE_DEVICES={mig} CUDA_MPS_PIPE_DIRECTORY={os.path.join(mps_pipe_dir, mig)} nvidia-cuda-mps-control'
        run_command(command)

def get_time(timer):
    return str(datetime.datetime.now()-timer).split('.', 1)[0]

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('MIG')
    parser.add_argument('--MPS', type=str, default=None)
    parser.add_argument('--kill_mps', type=bool, default=False)
    #parser.add_argument('--skip_mps_init', type=bool, default=False)
    parser.add_argument('--refresh', type=float, default=0.2)
    #parser.add_argument('--usev', type=str, default='')
    parser.add_argument('--query', type=str, default='/home/mdl/WGA_tests/query_blocked/')
    parser.add_argument('--target', type=str, default='/home/mdl/WGA_tests/target_blocked/')
    parser.add_argument('--tmp_dir', type=str, default='/home/mdl/WGA_tests/tmp/')
    parser.add_argument('--output', type=str, default='/home/mdl/WGA_tests/results/MIG_results.maf')
    parser.add_argument('--format', type=str, default='maf')
    parser.add_argument('--mps_pipe_dir', type=str, default='/home/mdl/WGA_tests/tmp/')
    parser.add_argument('--num_threads', type=int, default=-1)
    parser.add_argument('--segment_size', type=int, default=0)
    parser.add_argument('--segalign_cmd', type=str, default='run_segalign_symlink')
    parser.add_argument('--opt_cmd', type=str, default='')
    parser.add_argument('--test', type=bool, default=False)
    parser.add_argument('--keep_partial', type=bool, default=False)
    parser.add_argument('--only_missing', type=bool, default=False)
    parser.add_argument('--skip_mps_control', type=bool, default=False)




    args = parser.parse_args()
    
    if not os.path.exists(args.tmp_dir):
        print(f'Making directory for tmp_dir: {args.tmp_dir}')
        os.mkdir(args.tmp_dir)
    
    if not os.path.exists(args.mps_pipe_dir):
        print(f'Making directory for mps_pipe_dir: {args.mps_pipe_dir}')
        os.mkdir(args.mps_pipe_dir)
    
    output_format = args.output.split('.')[-1]

    mig_list = args.MIG.split(',')
    if len(mig_list) != len(set(mig_list)):
        print("Duplicate MIG ids given")
        assert(False)
    no_MIG = False
    if 'GPU' in args.MIG:
        print('Using non-MIG GPU')
        no_MIG = True

    if args.num_threads == -1:
        args.num_threads = os.cpu_count()
    print(f"USING {args.num_threads} THREADS")

    use_MPS = args.MPS != None
    if use_MPS:
        print("USING MPS")
        num_MPS = [int(x) for x in args.MPS.split(',')] # how many processes per MPS node
        max_proc_dict = dict(zip(mig_list, num_MPS))

    nvmlInit()

    test = False


    if test:
        print("PROCESSES")
        print(get_processes())
        print()
        print("UUIDS")
        print(get_uuids())
        print()
        print("MIG_DICT")
        mig_dict = init_mig_dict(mig_list)
        print(mig_dict)
        if no_MIG:
            print('NO_MIG DICT')
            mig_dict = {}
            mig_dict['0', 'N/A', 'N/A'] = mig_list[0]
            print(mig_dict)

        print("FREE MIG v1")
        print(get_free_mig(mig_dict))
        print("FREE MIG v2")
        print(get_free_mig_v2(mig_dict, no_MIG))
        input('...')
        while(True):
            print('enter loop')
            time.sleep(0.5)
            print(get_free_mig_v2(mig_dict, no_MIG))

    timer = datetime.datetime.now()
    cpu_per_segalign = 2
    gpu_per_segalign = 1
    #total_cpus = 128
    non_mig_gpu_id = 0 # TODO make this non hardcoded

    query_dir = args.query
    target_dir = args.target
    tmp_dir = os.path.abspath(args.tmp_dir)
    output_file = args.output

    mps_pipe_dir = os.path.abspath(args.mps_pipe_dir)

    uid_prefix = 'UID_'
    
    segment_size = args.segment_size

    segment_command = ''
    if args.segment_size != 0:
        segment_command = f'--segment_size {segment_size}'

    mps_timer = datetime.datetime.now()
    
    if args.kill_mps:
        destroy_mps(mig_list, mps_pipe_dir)
        exit(0)
        
    
    if use_MPS and not bool(int(args.skip_mps_control)):
        try:
            destroy_mps(mig_list, mps_pipe_dir)
        except:
            pass
        init_mps(mig_list, mps_pipe_dir)
    print(f"MPS init time: {get_time(mps_timer)}")

    _2bit_extension = '.2bit' #TODO make input parameter
    query_block_file_names = sorted([file for file in os.listdir(query_dir) if _2bit_extension not in file])
    target_block_file_names = sorted([file for file in os.listdir(target_dir) if _2bit_extension not in file])

    process_list = Process_List()#[]

    python_log = ""
    python_log += f"Starting Time: {datetime.datetime.now()}\n"

    normal_completion = False # used to check if program stopped due to error or termination
    
    resub_mem_fails = True

    try: # make sure to run exit functions on exit, even if errors
        if no_MIG:
            mig_dict = {}
            mig_dict[str(non_mig_gpu_id), 'N/A', 'N/A'] = mig_list[0]
        else:
            mig_dict = init_mig_dict(mig_list)

        print(mig_dict)

        part = 1
        refresh_time = args.refresh
        if not use_MPS:
            # TODO this if block is unneccessary, remove
            free_mig_list = get_free_mig_v2(mig_dict, no_MIG)
            for q in query_block_file_names:
                for t in target_block_file_names:
                    part += 1
                    #if part >= 10:
                    #    break
                    print(f"running part {part}: {q} and {t}")
                    while(len(free_mig_list) == 0):
                        time.sleep(refresh_time)
                        free_mig_list = get_free_mig_v2(mig_dict, no_MIG)

                    #command = f'CUDA_VISIBLE_DEVICES={free_mig_list[0]} run_segalign {os.path.join(target_dir, t)} {os.path.join(query_dir, q)} --debug --output {os.path.join(output_dir, f"{part}.maf")} --num_gpu {gpu_per_segalign} --num_threads {cpu_per_segalign} --num_lastz_threads {total_cpus}'
                    # TODO add CPU count parameter for lastz and segalign
                    command = f'CUDA_VISIBLE_DEVICES={free_mig_list[0]} run_segalign {os.path.join(target_dir, t)} {os.path.join(query_dir, q)} --debug --output={os.path.join(tmp_dir, f"part_{part}.{output_format}")} --num_gpu {gpu_per_segalign}'

                    # run process in non blocking way
                    process = NamedPopen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True, shell = True, name=f"segalign:{free_mig_list[0]}")
                    print(f"running process with pid={process.pid} and mig_uuid={free_mig_list[0]}. part {part}: {q} and {t}")
                    process_list.append(process)
                    free_mig_list.pop(0)
                    #assert(len(get_processes_v2()) <= len(mig_list)) # TODO why does this fail sometimes? => SegAlign spawns 2 CUDA processes per execution?
        else: # using MPS
            gpu_queue = GPU_queue(mig_list, tmp_dir, uid_prefix, num_MPS)
            print(f"GPU QUEUE: {gpu_queue.queue}")
            print("USING MPS SUBMISSION V3")

            mig_process_dict = gpu_queue.get_queue()
            free_mig_list = gpu_queue.get_free_device_list()

            max_processes = sum(max_proc_dict.values())
            print(f'max processes = {max_processes}')
            print(mig_process_dict)
            print('=====')
            pairs = [] # list of tasks
            for q in query_block_file_names:
                for t in target_block_file_names:
                    pairs.append((q,t))
            total_pairs = len(pairs)
            uid_pair_map = {} # used when resubmitting mem fails. Need to know file pair for given UID
            while len(pairs) > 0:
                #print(f"running part {part}: {q} and {t}")
                entered_while = False # used for printing to console
                while(len(free_mig_list) == 0):
                    entered_while = True
                    time.sleep(refresh_time)
                    gpu_queue.check_completion() # checks removes processes from gpu_queue whose UID files have been written

                    # some processes fail, need to deal with these # TODO prevent failures
                    running_uids = gpu_queue.get_running_uids()
                    failed_uids = process_list.check_uid_completion(running_uids) # uids that have their corresponding process complete, but not removed (possibly due to errors).
                                                                                  # Note: process removal done by checking if UID file is written
                    if resub_mem_fails:
                        resub_uids = process_list.check_uid_memfail(failed_uids)
                        for resub_uid in resub_uids:
                            print(f"{resub_uid} RESUBMITTED")
                            python_log += f"{resub_uid} RESUBMITTED\n"
                            resub_pair = uid_pair_map[resub_uid]
                            pairs.append(resub_pair)
                    #failed_uids = process_list.get_fails_and_check_completion()
                    if len(failed_uids) > 0:
                        gpu_queue.remove_uids(failed_uids)
                        for uid in failed_uids:
                            python_log += f"FAILED UID {uid}\n"
                            print(f"FAILED UID {uid}")

                    mig_process_dict = gpu_queue.get_queue()
                    free_mig_list = gpu_queue.get_free_device_list()
                    #print(f"Free devices: {free_mig_list}, device_queue: {gpu_queue.get_queue()}", end='\r')
                    #print(f"\33[2KFree devices: {free_mig_list}", end='\r')
                    print('\33[2K', end='') # clear line
                    
                    # print GPU queue information
                    for enum_i, i in enumerate(gpu_queue.get_queue().values()): # get single MIG slice queue
                        for enum_j,j in enumerate(i): # get
                            if enum_j == len(i)-1:
                                print(f"{j[4:]}", end='')
                            else:
                                print(f"{j[4:]} ", end='')
                        if enum_i == len(gpu_queue.get_queue())-1:
                            pass
                        else:
                            print(f"|", end='')
                    print(f"\r", end='')
                if entered_while:
                    print('')
                    
                # get MIG device
                mig_device, proc_ctr = free_mig_list[0] #TODO choose most remaining free space left device
                free_mig_list.pop(0)
                assert(proc_ctr < max_proc_dict[mig_device])
                while proc_ctr < max_proc_dict[mig_device]:
                    proc_ctr += 1
                    if len(pairs) == 0:
                        break
                    q,t = pairs[0]
                    pairs.pop(0)
                    uid_name = uid_prefix+str(part)
                    
                    out = os.path.join(tmp_dir, f"part_{part}.{output_format}")
                    command = f'CUDA_VISIBLE_DEVICES={mig_device} {args.segalign_cmd} {args.opt_cmd} {os.path.join(target_dir, t)} {os.path.join(query_dir, q)} --debug --output={out} --format={args.format} --num_gpu {gpu_per_segalign} --num_threads {args.num_threads} --uid {os.path.abspath(os.path.join(tmp_dir,uid_name))} {segment_command}'
                    command = f'CUDA_MPS_PIPE_DIRECTORY={os.path.join(mps_pipe_dir, mig_device)} ' + command
                    
                    if bool(args.only_missing):
                        if os.path.isfile(out):
                            part += 1
                            print(f"SKIPPING process, part {part}: {t} and {q}. Elapsed Time: {get_time(timer)}") 
                            python_log += 'SKIPPED: ' + command + '\n'
                            continue
                    # run process in non blocking way
                    process = NamedPopen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True, shell = True, name=uid_name)
                    process_list.append(process)
                    gpu_queue.submit(uid_name, mig_device)
                    running = sum([len(i) for i in gpu_queue.get_queue().values()])
                    est_runtime = str((datetime.datetime.now()-timer)*(total_pairs/ max(total_pairs-len(pairs)-running, 1)  )).split('.', 1)[0]
                    print(f"running process with pid={process.pid}, uid={uid_name} and mig_uuid={mig_device}. part {part} /{total_pairs}: {t} and {q}. Elapsed Time: {get_time(timer)}, estimated runtime: {est_runtime} [len(pairs) {len(pairs)}, running {running}]")
                    print(command)
                    #process_list.append(process)
                    mig_process_dict = gpu_queue.get_queue()
                    uid_pair_map[uid_name] = (q,t)
                    part += 1
                    python_log += command + '\n'

            # wait for all processes to complete
            output_stdout, output_stderr = process_list.get_output() #TODO waits here
            # TODO make sure that number of part files are expected
            print()
            print(f"Finished GPU section. time {get_time(timer)}")
            python_log += f"Finished GPU section. time {get_time(timer)}\n"


        # check if missing parts
        output_file_list = [i for i in os.listdir(tmp_dir) if i.startswith('part_') and i.endswith(f'.{output_format}')] 
        expected_outputs = len(query_block_file_names) * len(target_block_file_names)
        if len(output_file_list) != expected_outputs:
            print(f"Missing {expected_outputs-len(output_file_list)} output parts: ")
            set_output_file_list = set(output_file_list)
            
            for i in range(1, expected_outputs+1):
                part_name = f"part_{i}.{output_format}"
                if part_name not in set_output_file_list:
                    print(f"{i}, ", end='')
            print()
            print("Rerun with --only_missing")
        
        normal_completion = True
    finally:

        if not normal_completion:
            output_stdout, output_stderr = process_list.get_output(terminate=True)
        stdout_log_file = process_list.stdout
        stderr_log_file = process_list.stderr
        # log outputs to tmp directory
        stdout_log_file = 'o.txt'#f'{tmp_dir}o.txt'
        write_mode = 'w'
        if args.only_missing:
            write_mode = 'a'
        with open(stdout_log_file, write_mode) as file:
            print(f"Writing stdout to {stdout_log_file}")
            file.write(output_stdout)
        stderr_log_file = 'e.txt'#f'{tmp_dir}e.txt'
        with open(stderr_log_file, write_mode) as file:
            print(f"Writing stderr to {stderr_log_file}")
            file.write(output_stderr)


        print(f"total time {get_time(timer)}")
        python_log += f"total time {get_time(timer)}\n"
        combine_results(tmp_dir, output_file, part_pattern=f'part_*.{output_format}', remove=(not bool(args.keep_partial)))
        print(f"result combination finished at {get_time(timer)}")
        python_log += f"result combination finished at {get_time(timer)}\n"
        python_log += f"Ending Time: {datetime.datetime.now()}\n"
        # finished execution
        python_log_file = 'p.txt'#f'{tmp_dir}p.txt'
        with open(python_log_file, write_mode) as file:
            print(f"Writing executed commands to {python_log_file}")
            file.write(python_log)
        #print(output_stderr)
        #print(output_stdout)
        
        if not bool(not bool(args.keep_partial)):
            remove_uids(tmp_dir, uid_prefix)
        if use_MPS and not bool(int(args.skip_mps_control)):
            destroy_mps(mig_list, mps_pipe_dir)

        nvmlShutdown()
