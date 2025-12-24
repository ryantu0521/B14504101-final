import numpy as np
import matplotlib.pyplot as plt

COMMON_ACIDS = {
    "1": {"name": "鹽酸 (HCl - 強酸)", "eng_name": "HCl", "pKas": [-7]},
    "2": {"name": "硝酸 (HNO3 - 強酸)", "eng_name": "HNO3", "pKas": [-1.4]},
    "3": {"name": "硫酸 (H2SO4 - 強酸+弱酸)", "eng_name": "H2SO4", "pKas": [-3, 1.99]},
    "4": {"name": "醋酸 (CH3COOH - 弱酸)", "eng_name": "Acetic Acid", "pKas": [4.76]},
    "5": {"name": "草酸 (H2C2O4 - 二質子酸)", "eng_name": "Oxalic Acid", "pKas": [1.25, 4.27]},
    "6": {"name": "碳酸 (H2CO3 - 二質子酸)", "eng_name": "Carbonic Acid", "pKas": [6.35, 10.33]},
    "7": {"name": "磷酸 (H3PO4 - 三質子酸)", "eng_name": "Phosphoric Acid", "pKas": [2.12, 7.21, 12.32]},
    "8": {"name": "氫氟酸 (HF - 弱酸)", "eng_name": "HF", "pKas": [3.17]},
    "0": {"name": ">> 自訂酸 (Custom) <<", "eng_name": "Custom Acid", "pKas": []}
}

COMMON_BASES = {
    "1": {"name": "氫氧化鈉 (NaOH - 強鹼)", "eng_name": "NaOH", "pKbs": []}, # 空列表代表強鹼
    "2": {"name": "氫氧化鉀 (KOH - 強鹼)", "eng_name": "KOH", "pKbs": []},
    "3": {"name": "氨水 (NH3 - 弱鹼)", "eng_name": "Ammonia", "pKbs": [4.75]},
    "4": {"name": "碳酸鈉 (Na2CO3 - 二元鹼)", "eng_name": "Sodium Carbonate", "pKbs": [3.67, 7.65]},
    "5": {"name": "乙二胺 (en - 二元鹼)", "eng_name": "Ethylenediamine", "pKbs": [3.29, 6.44]},
    "6": {"name": "三乙胺 (Triethylamine - 弱鹼)", "eng_name": "Triethylamine", "pKbs": [3.25]},
    "7": {"name": "吡啶 (Pyridine - 弱鹼)", "eng_name": "Pyridine", "pKbs": [8.77]},
    "0": {"name": ">> 自訂鹼 (Custom) <<", "eng_name": "Custom Base", "pKbs": []}
}

INDICATORS = [
    {"name": "Methyl Orange (甲基橙)", "low": 3.1, "high": 4.4, "color": "orange"},
    {"name": "Methyl Red (甲基紅)", "low": 4.4, "high": 6.2, "color": "red"},
    {"name": "Bromothymol Blue (溴瑞香草酚藍)", "low": 6.0, "high": 7.6, "color": "green"},
    {"name": "Phenolphthalein (酚酞)", "low": 8.2, "high": 10.0, "color": "magenta"},
    {"name": "Thymolphthalein (百里酚酞)", "low": 9.3, "high": 10.5, "color": "blue"},
    {"name": "Alizarin Yellow (茜素黃)", "low": 10.1, "high": 12.0, "color": "yellow"}
]

def calculate_acid_mean_charge(h, pKas):
    """計算酸的平均負電荷"""
    if not pKas: return 1.0 
    Kas = [10**(-p) for p in pKas]
    n = len(Kas)
    denom = 0.0
    h_powers = [h**(n - i) for i in range(n + 1)]
    k_cum = [1.0]
    temp_k = 1.0
    for k in Kas:
        temp_k *= k
        k_cum.append(temp_k)
    for i in range(n + 1):
        denom += h_powers[i] * k_cum[i]
    if denom == 0: return 0.0
    mean_charge = 0.0
    for i in range(1, n + 1):
        term = k_cum[i] * (h ** (n - i))
        alpha = term / denom
        mean_charge += i * alpha
    return mean_charge

def calculate_base_mean_charge(target_ph, pKbs):
    if not pKbs: return 1.0 
    Kw = 1e-14
    h = 10**(-target_ph)
    if h == 0: h = 1e-20
    oh = Kw / h
    Kbs = [10**(-p) for p in pKbs]
    n = len(Kbs)
    denom = 0.0
    oh_powers = [oh**(n - i) for i in range(n + 1)]
    k_cum = [1.0]
    temp_k = 1.0
    for k in Kbs:
        temp_k *= k
        k_cum.append(temp_k)
    for i in range(n + 1):
        denom += oh_powers[i] * k_cum[i]
    if denom == 0: return 0.0
    mean_charge = 0.0
    for i in range(1, n + 1):
        term = k_cum[i] * (oh ** (n - i))
        alpha = term / denom
        mean_charge += i * alpha
    return mean_charge

def solve_general_volume(target_ph, mode, acid_data, acid_conc, base_data, base_conc, beaker_vol):
    Kw = 1e-14
    h = 10 ** (-target_ph)
    if h == 0: h = 1e-20
    oh = Kw / h
    
    Q_acid = calculate_acid_mean_charge(h, acid_data['pKas'])
    Q_base = calculate_base_mean_charge(target_ph, base_data['pKbs']) 

    water_term = h - oh
    
    if mode == 'A': 
        numerator = beaker_vol * (Q_acid * acid_conc - water_term)
        denominator = Q_base * base_conc + water_term
    else: 
        numerator = beaker_vol * (Q_base * base_conc + water_term)
        denominator = Q_acid * acid_conc - water_term
    
    if abs(denominator) < 1e-15: return -1.0
    return numerator / denominator

def get_user_choice(options_dict, prompt_msg):
    print(f"\n--- {prompt_msg} ---")
    keys = sorted(options_dict.keys(), key=lambda x: int(x))
    for key in keys:
        val = options_dict[key]
        print(f"  {key}. {val['name']}")
    while True:
        choice = input(f"請選擇 (輸入代號): ").strip()
        if choice in options_dict: return choice
        print("無效的代號，請重新輸入。")

def get_float_input(prompt):
    """強制輸入有效浮點數"""
    while True:
        user_input = input(f"{prompt}: ").strip()
        if not user_input:
            print("錯誤：請勿留空，請輸入數值。")
            continue
        try:
            value = float(user_input)
            if value <= 0:
                print("錯誤：數值必須大於 0。")
                continue
            return value
        except ValueError:
            print("錯誤：輸入格式不正確，請輸入數字。")

def check_weak_vs_weak(acid_data, base_data):
    """
    檢查是否為 弱酸 vs 弱鹼。
    定義：
      弱酸: 第一個 pKa > 0 (強酸通常是負的)
      弱鹼: 有 pKb 值 (強鹼通常為空列表或極小)
    """
    is_weak_acid = False
    if acid_data['pKas'] and acid_data['pKas'][0] > 0:
        is_weak_acid = True
    
    is_weak_base = False
    if base_data['pKbs'] and len(base_data['pKbs']) > 0:
        is_weak_base = True
        
    if is_weak_acid and is_weak_base:
        return True
    return False

def main():
    print("      酸鹼滴定模擬器      ")
    while True:
        mode_input = input("\n請選擇滴定方向 (A: 酸在燒杯 / B: 鹼在燒杯): ").strip().upper()
        if mode_input not in ['A', 'B']:
            print("無效輸入，請輸入 A 或 B。")
            continue
        
        mode = mode_input
        
        if mode == 'A':
            acid_choice = get_user_choice(COMMON_ACIDS, "選擇酸 (Analyte - 燒杯)")
            base_choice = get_user_choice(COMMON_BASES, "選擇鹼 (Titrant - 滴定管)")
        else:
            base_choice = get_user_choice(COMMON_BASES, "選擇鹼 (Analyte - 燒杯)")
            acid_choice = get_user_choice(COMMON_ACIDS, "選擇酸 (Titrant - 滴定管)")

        temp_acid = COMMON_ACIDS[acid_choice].copy()
        temp_base = COMMON_BASES[base_choice].copy()

        if acid_choice == '0':
            print("\n[自訂酸設定]")
            temp_acid['eng_name'] = input("請輸入酸的英文簡稱: ")
            while True:
                try:
                    temp_acid['pKas'] = sorted([float(x) for x in input("請輸入 pKas (逗號分隔): ").split(',')])
                    break
                except: print("格式錯誤。")
        
        if base_choice == '0':
            print("\n[自訂鹼設定]")
            temp_base['eng_name'] = input("請輸入鹼的英文簡稱: ")
            yn = input("是強鹼嗎? (y/n): ").lower()
            if yn == 'y': temp_base['pKbs'] = []
            else:
                while True:
                    try:
                        temp_base['pKbs'] = sorted([float(x) for x in input("請輸入 pKbs (逗號分隔): ").split(',')])
                        break
                    except: print("格式錯誤。")

        if check_weak_vs_weak(temp_acid, temp_base):
            print("\n" + "="*50)
            print("錯誤：本系統不支援「弱酸」與「弱鹼」互相滴定！")
            print("原因：此組合無明顯當量點，無法進行有效滴定分析。")
            print("請重新選擇含有至少一個強電解質的組合。")
            print("="*50)
            continue
        else:
            acid_data = temp_acid
            base_data = temp_base
            break

    analyte_data = acid_data if mode == 'A' else base_data
    titrant_data = base_data if mode == 'A' else acid_data
    
    print("\n--- 設定實驗參數 (請輸入數值) ---")
    conc_analyte = get_float_input(f"請輸入 燒杯[{analyte_data['eng_name']}] 的濃度 (M)")
    vol_analyte  = get_float_input(f"請輸入 燒杯[{analyte_data['eng_name']}] 的體積 (mL)")
    conc_titrant = get_float_input(f"請輸入 滴定管[{titrant_data['eng_name']}] 的濃度 (M)")

    n_acid_protons = len(acid_data['pKas']) if acid_data['pKas'] else 1
    n_base_protons = len(base_data['pKbs']) if base_data['pKbs'] else 1
    
    if mode == 'A':
        acid_conc, base_conc, beaker_vol = conc_analyte, conc_titrant, vol_analyte
        total_protons = n_acid_protons
        equiv_vol = (conc_analyte * vol_analyte * total_protons) / conc_titrant
    else:
        acid_conc, base_conc, beaker_vol = conc_titrant, conc_analyte, vol_analyte
        total_protons = n_base_protons
        equiv_vol = (conc_analyte * vol_analyte * total_protons) / conc_titrant

    analyte_name = analyte_data['eng_name']
    titrant_name = titrant_data['eng_name']

    ph_range = np.linspace(0, 14, 2000)
    v_results = []
    ph_results = []
    max_plot_vol = equiv_vol * 3.5

    for ph in ph_range:
        v = solve_general_volume(ph, mode, acid_data, acid_conc, base_data, base_conc, beaker_vol)
        if 0 <= v <= max_plot_vol:
            v_results.append(v)
            ph_results.append(ph)

    v_arr = np.array(v_results)
    ph_arr = np.array(ph_results)
    
    eq_points = []
    
    if mode == 'A':
        num_pts = len(acid_data['pKas']) if acid_data['pKas'] else 1
        for i in range(1, num_pts + 1):
            v_eq = (conc_analyte * vol_analyte * i) / conc_titrant
            if v_eq <= max_plot_vol:
                idx = (np.abs(v_arr - v_eq)).argmin()
                eq_points.append((v_eq, ph_arr[idx], i))
    else:
        num_pts = len(base_data['pKbs']) if base_data['pKbs'] else 1
        for i in range(1, num_pts + 1):
            v_eq = (conc_analyte * vol_analyte * i) / conc_titrant
            if v_eq <= max_plot_vol:
                idx = (np.abs(v_arr - v_eq)).argmin()
                eq_points.append((v_eq, ph_arr[idx], i))

    best_indicator = None
    if eq_points:
        last_ph = eq_points[-1][1]
        candidates = [ind for ind in INDICATORS if (ind['low']-0.5) <= last_ph <= (ind['high']+0.5)]
        if candidates:
            best_indicator = min(candidates, key=lambda x: abs((x['low']+x['high'])/2 - last_ph))
            print(f"\n推薦指示劑: {best_indicator['name']} (變色範圍: {best_indicator['low']}-{best_indicator['high']})")
        else:
            print("\n未找到完美匹配的常見指示劑")

    plt.figure(figsize=(10, 7))
    color = 'navy'
    
    if best_indicator:
        plt.axhspan(best_indicator['low'], best_indicator['high'], color=best_indicator['color'], alpha=0.15)
        indicator_eng_name = best_indicator['name'].split('(')[0].strip()
        plt.text(max_plot_vol*0.02, best_indicator['high']+0.2, 
                 f"Indicator: {indicator_eng_name}", 
                 color='black', fontsize=9)

    plt.plot(v_arr, ph_arr, color=color, linewidth=2.5, label='pH Curve')

    for v_eq, ph_eq, idx in eq_points:
        plt.scatter(v_eq, ph_eq, color='red', s=60, zorder=10)
        
        is_polyprotic = len(eq_points) > 1
        if is_polyprotic:
            label_text = f"Eq Point {idx}\n({v_eq:.1f} mL)"
        else:
            label_text = f"Eq Point\n({v_eq:.1f} mL)"
            
        plt.text(v_eq, ph_eq + 0.8, label_text, color='red', ha='center', fontweight='bold', fontsize=10)
        plt.axvline(x=v_eq, color='gray', linestyle=':', alpha=0.5)

    plt.title(f'Titration: {analyte_name} vs {titrant_name}', fontsize=16)
    plt.xlabel(f'Volume of {titrant_name} Added (mL)', fontsize=14)
    plt.ylabel('pH', fontsize=14)
    plt.legend(loc='best')
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.ylim(0, 14)
    plt.xlim(0, max_plot_vol)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
