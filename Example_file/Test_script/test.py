'''
Descripttion: 
Author: Ne0tea
version: 
Date: 2024-09-27 00:00:05
LastEditors: Ne0tea
LastEditTime: 2024-09-27 16:27:03
'''
def split_string_every_4_chars(s):
    ### Should be aligned to collapse format
    cur_list=s.split("|")
    cur_list = [x for x in cur_list if x]
    # cur_list=[s[i:i+4] for i in range(0, len(s), 4)]
    return len(list(set(cur_list)))

def any_key_less(dict1, dict2):
    for key, value in dict2.items():
        if key in dict1 :
            dict1_vlaue=split_string_every_4_chars(dict1[key])
            if value > dict1_vlaue:
                return True
        else:
            return True
    return False

cur={'BB': 'Bamp|DlatC|Phedu|Bamp|Olat|Bamp|DlatB|DlatA|Bamp|Phedu|Olat|Phedu', 'OO': 'Lper|Obra|Opun|Zpal|Obar|Oglu|Omer|Ogla|Zlat|Oruf|Osat|Zlat|Zpal', 'PP': 'Aspl|Bdis|Bsta|Aspl|Bdis|Bsta', 'PAC': 'Ehap|EcruA|EoryA|EcruC|EcolD|EcruB|MsinB|MsinA|Sbic|Zmay|EcolD|Doli|DexiB|EcruA|Sita|Ehap|Asem|EoryB|EoryA|Phal|EcruC|Svir|PvirN|EcolF|EcolE', 'MAD': 'EtefB|CsonA|Enin|CsonB|Zjap|Ecur|Otho|EtefA', 'OUT': ''}
ss={'OUT':'Plat'}
temp={'BB': 2, 'OO': 3, 'PP': 4, 'PAC': 6, 'MAD': 3, 'OUT': 1}
xxx=(1,2,8,2,3,3,4)
print(len(xxx))
print(any_key_less(cur,temp))