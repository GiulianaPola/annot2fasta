#!/usr/bin/env python3
# -*- coding: utf-8 -*-
entry_id=None
try:
    import traceback
    import os
    import argparse
    import datetime
    import sys
    import re
    
    
    version="1.1.0"
    print('annot2fasta {}'.format(version))
    
    start_time = datetime.datetime.now()
    call=os.path.abspath(os.getcwd())
    
    annot_files=dict()
    error=[]
    
    parser = argparse.ArgumentParser(
      prog='annot2fasta',
      description=
      'Create FASTA files of nucleotide sequences from annotation files.',
      epilog='https://github.com/GiulianaPola/annot2fasta',add_help=False)
    parser.add_argument('-i',help="Annotation file or folder")
    parser.add_argument('-o',help='Name of the output folder',default="fasta_files")
    parser.add_argument('-h', '-help', action='store_true')
    parser.add_argument('-v', action='store_true',help='Version')
    #parser.add_argument('-config', help='Configuration file')
    args=parser.parse_args()
    
    def printing_help():
        print('(c) 2023. Giuliana Pola & Arthur Gruber\n')
        print('For more information access: ')
        print('Usage: annot2fasta -i <path|file>')
        print('Create FASTA files of nucleotide sequences from annotation files.')
        print('\nOptional arguments:')
        print('  -h, --help             Show this help message and exit')
        print('  -i <path|file>        Annotation files or folder (required)')
        print('  -o string       Name of the output folder (default: "fasta_files")')
        print('  -v                    Show version information')
    
    def rename(name):
     i=0
     path = os.path.dirname(name)
     name = os.path.basename(name)
     newname = os.path.join(path, name)
     while os.path.exists(newname):
         i += 1
         newname = os.path.join(path, "{}_{}".format(name, i))
     return newname
    
    def validate(args):
     args.i=os.path.abspath(args.i)
     print("\nargs.i={}\n".format(args.i))
     if os.path.isdir(args.i):
       if not args.i[-1]=="/":
         args.i=args.i+"/"
       if not os.path.exists(args.i):
           print("ERROR: Annotation files' folder '{}' does not exist!".format(args.i))
           exit()
       elif not os.listdir(args.i):
         print("ERROR: Annotation files' folder '{}' is empty!".format(args.i))
         exit()
       elif os.path.isdir(args.i):
         for root, dirs, files in os.walk(args.i):
           for file in files:
               if file.endswith('.tab'):
                   file=os.path.join(root, file)
                   if os.path.split(file)[0]+"/" not in annot_files:
                     annot_files[os.path.split(file)[0]+"/"]=[]
                   annot_files[os.path.split(file)[0]+"/"].append(os.path.split(file)[-1])
         if annot_files==[]:
           print("ERROR: Input folder '{}' does not contain annotation files (.tab)!".format(args.i))
           exit()
     elif os.path.isfile(args.i):
       if not os.path.exists(args.i):
           print("ERROR: Annotation file '{}' does not exist!".format(args.i))
           exit()
       else:
         annot_files[os.path.split(args.i)[0]+"/"]=[os.path.split(args.i)[-1]]
       
       
     try:
       if os.path.exists(args.o):
         args.o=rename(args.o)
       os.mkdir(args.o)
       args.o=os.path.abspath(args.o)
     except Exception as e:
       print(e)
       print("ERROR: Output folder is not writable!")
       exit()
     return args.i
          
    
    def odd_list(coord, start, end):
        log=open(os.path.join(args.o, 'file.log'), 'a')
        dif = dict()
        for c in coord:
            a=int(c[0])
            b=int(c[-1])
            temp1 = [s for s in start if s != a]
            temp2 = [e for e in end if e != b]
    
            dif_s = sum(abs(s - a) for s in temp1) if temp1 else 0
            dif_e = sum(abs(e - b) for e in temp2) if temp2 else 0
    
            key = ((dif_s / len(temp1) if temp1 else 0) + (dif_e / len(temp2) if temp2 else 0)) / 2
    
            if key not in dif:
                dif[key] = []
            dif[key].append([a, b])
        max_key = max(dif.keys())
        return dif[max_key]
    
    def split_list(coord, med):
        size=int(len(coord)/2)
        result = []
        i = 1
        temp = []
    
        for j in range(len(med)):
            if med[j] <= ((max(med) - min(med)) / size * i) + min(med):
                temp.append(coord[j])
            else:
                result.append(temp)
                temp = [coord[j]]
                i = i + 1
            result.append(temp)
    
        return result
    
    def create_fasta_files(annot_files, fasta_folder):
        for annot_path in annot_files:
            if os.path.isdir(args.i) and not annot_path==args.i:
              out_path=os.path.join(fasta_folder,annot_path.replace(args.i,""))
              log.write("\nFolder: {}\n".format(annot_path.replace(args.i,"")))
              if not os.path.exists(out_path):
                os.makedirs(out_path)
            else:
              out_path=fasta_folder
            for annot in annot_files[annot_path]:
              extract_sequences(os.path.join(annot_path,annot),out_path)
        if error==[]:
          log.write("\nFASTA files created successfully!\n")
          print("FASTA files created successfully!")
        else:
          if len(error)==1:
            log.write("\nERROR: {} FASTA file was not created successfully: {}\n".format(len(error),error))
            print("ERROR: {} FASTA file was not created successfully: {}".format(len(error),error))
          else:
            log.write("\nERROR: {} FASTA files was not created successfully: {}\n".format(len(error),error))
            print("ERROR: {} FASTA files was not created successfully: {}".format(len(error),error))
    
    def remove_duplicates(coord):
        set_lists = set(map(tuple, coord))
        unique_lists = [list(t) for t in set_lists]
        return sorted(unique_lists)
    
    def extract_sequences(annot, out_path):
        coordinates = []
        sequence = ''
        entry_id = os.path.split(annot)[-1]
        if "_all_results.tab" in entry_id:
            entry_id = entry_id.replace("_all_results.tab", "")
    
        try:
            annot=os.path.abspath(annot)
            out_path=os.path.abspath(out_path)
            log.write("\nannotation file: {}\n".format(annot))
            log.write("fasta folder: {}\n".format(out_path))
            with open(annot, "r") as annot_file:
                repeat_region = []
                element = []
                lines = annot_file.readlines()
                j=0
                for i in range(len(lines)):
                    search=''
                    line = lines[i].strip()
    
                    if '/label=element' in line:
                        search=''
                        for o in range(i, j, -1):
                          temp = lines[o]
                          if "/" not in temp and "=" not in temp:
                            search=temp
                            break
                        if not search=='':
                          search=search.strip().split()[-1]
                          if "join(" in search:
                            search=search.replace("join(","").replace(")","")
                          if ".." in search and len(search.split("..")) == 2:
                              coords=search.split("..")
                              if not [int(coords[0]), int(coords[1])] in element:
                                element.append([int(coords[0]), int(coords[1])])
                                log.write("{} {}: {}\n".format(entry_id, line.split("/label=")[-1].strip(), search))
                          else:
                              log.write("'{}' WARNING: Skipping {} invalid coordinate format: '{}'\n".format(entry_id, line.split("/label=")[-1].strip(),search))
                    if '/note="direct repeat"' in line or '/note="inverted repeat"' in line or '/rpt_type=inverted' in line or '/standard_name="TDR"' in line or '/standard_name="TIR"' in line:
                        search=''
                        for o in range(i, j, -1):
                          temp = lines[o]
                          if "/" not in temp and "=" not in temp:
                            search=temp
                            break
                        if not search=='':
                          search=search.strip().split()[-1]
                          if "join(" in search:
                            search=search.replace("join(","").replace(")","")
                          if "," in search and len(search.split(",")) == 2:
                              coords = search.split(",")
                              start = int(coords[0].split("..")[0])
                              end = int(coords[1].split("..")[1])
                              if not [start, end] in repeat_region:
                                repeat_region.append([start, end])
                                log.write("{} {}: {}\n".format(entry_id, line.split("=")[-1].strip(), search))
                          else:
                              log.write("'{}' WARNING: Skipping {} invalid coordinate format: '{}'\n".format(entry_id,line.split("=")[-1].strip(),search))
                    elif line.startswith('//'):
                        break
                    elif line.startswith('SQ') or line.startswith('FT'):
                        continue
                    else:
                        line = re.sub(r'\d+', '', line)
                        line = re.sub(r'\s', '', line)
                        sequence += line.upper()
    
            if not repeat_region == []:
                repeat_region = remove_duplicates(repeat_region)
                if not element == []:
                  within=True
                  for elm in element:
                    for rr in repeat_region:
                      if elm[0]<rr[0] or elm[-1]>rr[-1]:
                        within=False
                  if within==False:
                    log.write("{} WARNING: The previously annotated element is not within the repeat region!\n".format(entry_id))
                if len(repeat_region) > 2:
                    log.write("{} WARNING: '{}' has {} repeat regions!\n".format(entry_id, entry_id, len(repeat_region)))
                med = []
                dif = []
                result = []
                size = int(len(repeat_region) / 2)
                start = []
                end = []
                for rr in repeat_region:
                    start.append(int(rr[0]))
                    end.append(int(rr[-1]))
                    med.append((int(rr[0]) + int(rr[-1])) / 2)
                if len(repeat_region) % 2 == 1 and len(repeat_region)>2:
                    dif = odd_list(sorted(repeat_region), start, end)
                    if len(dif) == 1:
                        log.write("{} WARNING: Removing the outlier repeat_region {} from the list...\n".format(entry_id, dif[0]))
                        repeat_region.remove(dif[0])
                if len(repeat_region) <= 3:
                    for rr in repeat_region:
                        start.append(int(rr[0]))
                        end.append(int(rr[-1]))
                    coordinates.append([sorted(start)[0], sorted(end)[-1]])
                else:
                    for elmt in split_list(sorted(repeat_region), med):
                        start = []
                        end = []
                        for e in elmt:
                            start.append(int(e[0]))
                            end.append(int(e[-1]))
                        coordinates.append([sorted(start)[0], sorted(end)[-1]])
            elif not element == []:
                element = remove_duplicates(element)
                start = []
                end = []
                for e in element:
                    start.append(int(e[0]))
                    end.append(int(e[-1]))
                coordinates.append([sorted(start)[0], sorted(end)[-1]])
    
            if coordinates == []:
                log.write("'{}' ERROR: The coordinates of the '{}' element could not be found!\n".format(entry_id))
                print("'{}' ERROR: The coordinates of the element could not be found!\n".format(entry_id))
                error.append(entry_id)
                try:
                  if not os.path.isdir(os.path.join(out_path,"with_element")):
                      os.makedirs(os.path.join(out_path,"with_element"))
                  with open(os.path.join(out_path,"{}_with_element.fasta".format(entry_id)), "w") as fasta_file:
                      fasta_file.write(">{} With Element\n".format(entry_id))
                      log.write(">{} With Element\n".format(entry_id))
                      seq = sequence
                      for i in range(0, len(seq), 70):
                          fasta_file.write(seq[i:i + 70] + "\n")
                except Exception as e:
                   log.write("'{}' ERROR: {}\n".format(entry_id,e))
                   error.append(entry_id)

            else:
                try:
                  fasta = []
                  coordinates = remove_duplicates(coordinates)
                  for start, end in coordinates:
                      fasta.append("{}-{}".format(start, end))
                  if not os.path.isdir(os.path.join(out_path,"element")):
                      os.makedirs(os.path.join(out_path,"element"))
                  k = 1
                  for start, end in sorted(coordinates):
                      if k==1:
                        fasta_file="{}_element.fasta".format(entry_id)
                      elif entry_id[-2:]=="_1":
                        new_entry_id="{}_{}".format(entry_id[:-2],str(k))
                        fasta_file="{}_element.fasta".format(new_entry_id)
                      with open(os.path.join(out_path,"element",fasta_file), "w") as fasta_handler:
                          fasta_handler.write(">{} Element_{} {}-{}\n".format(entry_id, k, start, end))
                          log.write(">{} Element_{} {}-{}\n".format(entry_id, k, start, end))
                          seq = sequence[start - 1:end]
                          for i in range(0, len(seq), 70):
                              fasta_handler.write(seq[i:i + 70] + "\n")
                          k += 1
                except Exception as e:
                   print("'{}' ERROR: {}\n".format(entry_id, str(e)))
                   exc_type, exc_value, exc_traceback = sys.exc_info()
                   formatted_lines = traceback.format_exc().splitlines()
                   for i, line in enumerate(formatted_lines):
                      print("{}: {}".format(i, line))
                   error.append(entry_id)
                
                try:
                  if not os.path.isdir(os.path.join(out_path,"with_element")):
                      os.makedirs(os.path.join(out_path,"with_element"))
                  with open(os.path.join(out_path,"with_element","{}_with_element.fasta".format(entry_id)), "w") as fasta_file:
                      fasta_file.write(">{} with element in {}\n".format(entry_id, ",".join(fasta)))
                      log.write(">{} with element in {}\n".format(entry_id, ",".join(fasta)))
                      seq = sequence
                      for i in range(0, len(seq), 70):
                          fasta_file.write(seq[i:i + 70] + "\n")
                except Exception as e:
                   print("'{}' ERROR: {}\n".format(entry_id, str(e)))
                   exc_type, exc_value, exc_traceback = sys.exc_info()
                   formatted_lines = traceback.format_exc().splitlines()
                   for i, line in enumerate(formatted_lines):
                      print("{}: {}".format(i, line))
                   error.append(entry_id)
                
                try:
                  if not os.path.isdir(os.path.join(out_path,"without_element")):
                      os.makedirs(os.path.join(out_path,"without_element"))
                  with open(os.path.join(out_path,"without_element","{}_without_element.fasta".format(entry_id)), "w") as fasta_file:
                      minimum = 0
                      coordinates.append([len(sequence), 0])
                      fasta = []
                      seq = ''
                      last_end=-1
                      last_start=-1
                      for start, end in coordinates:
                          if not last_end==-1:
                            if start<last_end:
                              coordinates.remove([start,end])
                              coordinates.remove([last_start,last_end])
                              coordinates.append([last_start,end])
                          last_start=start
                          last_end=end
                      for start, end in sorted(coordinates):
                          start += -1
                          end += -1
                          fasta.append("{}-{}".format(minimum + 1, start))
                          if end == -1:
                              seq = seq + sequence[minimum:]
                          else:
                              seq = seq + sequence[minimum:start]
                          minimum = end + 1
                      if not minimum == 0:
                          log.write("'{}' ERROR: Error extracting element from sequence\n".format(entry_id))
                      fasta_file.write(">{} Without Element {}\n".format(entry_id, ",".join(fasta)))
                      log.write(">{} Without Element {}\n".format(entry_id, ",".join(fasta)))
                      for i in range(0, len(seq), 70):
                          fasta_file.write(seq[i:i + 70] + "\n")
                except Exception as e:
                   print("'{}' ERROR: {}\n".format(entry_id, str(e)))
                   exc_type, exc_value, exc_traceback = sys.exc_info()
                   formatted_lines = traceback.format_exc().splitlines()
                   for i, line in enumerate(formatted_lines):
                      print("{}: {}".format(i, line))
                   error.append(entry_id)
    
        except Exception as e:
            error.append(entry_id)
            print("'{}' ERROR: {}".format(entry_id, e))
            exc_type, exc_value, exc_traceback = sys.exc_info()
            formatted_lines = traceback.format_exc().splitlines()
            for i, line in enumerate(formatted_lines):
                print("{}: {}".format(i, line))
    
    if len(sys.argv)==1:
        printing_help()
        exit()
    if "-h" in sys.argv or "--help" in sys.argv:
        printing_help()
        exit()
    if args.v:
        print(version)
        exit()
    elif not args.i:
        print("ERROR: Missing annotation files' folder (-i)!")
        exit()
    
    arg_input = validate(args)
    
    try:
      log=open(os.path.join(args.o, 'file.log'), 'w')
      try:
        log.write('annot2fasta v{}\n'.format(version))
        log.write('(c) 2023. Giuliana Pola & Arthur Gruber\n')
        log.write('For more information access: ')
        log.write('\nStart time: {}\n'.format(start_time.strftime("%d/%m/%Y, %H:%M:%S")))
        log.write('\nWorking directory: {}\n'.format(call))
        log.write('\nCommand line: {}\n'.format(' '.join(sys.argv)))
        user=""
        try:
          user=os.getlogin()
        except Exception as e:
          try:
            user=os.environ['LOGNAME']
          except Exception as e:
            try:
              user=os.environ['USER']
            except Exception as e:
              pass
            else:
              pass
          else:
            pass
        else:
          pass
        if not user=="":
          log.write('\nUser: {}\n'.format(user))
        log.write('\nParameters:\n')
        for arg in vars(args):
            value = getattr(args, arg)
            if value is not None and value is not False:
                log.write("{}={}\n".format(arg,value))
      except Exception as e:
        print("ERROR: Log file was not written!")
        exit()
    except Exception as e:
        print("ERROR: Log file was not created!")
        exit()
    log.close()
    log=open(os.path.join(args.o, 'file.log'), 'a')
    create_fasta_files(annot_files, args.o)
    execution=datetime.datetime.now() - start_time
    print("\nExecution time: {}".format(execution))
    log.write("\nExecution time: {}\n".format(execution))
    print("\nExecution time per file: {}".format(execution/len(str(list(annot_files.values())).replace("[","").replace("]","").replace("'","").replace(", ",",").split(","))))
    log.write("\nExecution time per file: {}\n".format(execution/len(str(list(annot_files.values())).replace("[","").replace("]","").replace("'","").replace(", ",",").split(","))))
    print("End")
    log.write("\nEnd\n")
    log.close()
except Exception as e:
    if not entry_id==None:
      error.append(entry_id)
      print("ERROR '{}': {}".format(entry_id, e))
    else:
      print("ERROR: {}".format(e))
    exc_type, exc_value, exc_traceback = sys.exc_info()
    formatted_lines = traceback.format_exc().splitlines()
    for i, line in enumerate(formatted_lines):
      print("{}: {}".format(i, line))                


    

