from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import matplotlib.colors as mcolors
from matplotlib.colors import hsv_to_rgb

# 互补碱基对
complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# 生成反向互补序列
def reverse_complement(seq):
    return ''.join([complement[base] for base in reversed(seq)])

# 生成更美观的颜色
def generate_colors(n):
    colors = []
    for i in range(n):
        hue = (i * 0.618) % 1.0  # 使用黄金比例分布色相
        saturation = 0.8 if i % 2 else 0.6  # 交替饱和度
        value = 0.9 if i % 3 else 0.8  # 交替亮度
        rgb = hsv_to_rgb((hue, saturation, value))
        rgba = (*rgb, 0.8)  # 提高透明度
        hex_color = mcolors.rgb2hex(rgb)
        rgba_color = f"rgba({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)}, 0.8)"
        colors.append((hex_color, rgba_color))
    return colors

def calculate_gc_content(seq1, seq2, start1, end1, start2, end2):
    """计算配对区域中的GC含量"""
    # 获取实际配对区域
    region1 = seq1[start1:end1+1]
    region2 = seq2[start2:end2+1]
    
    # 统计所有G和C碱基
    gc_count = region1.count('G') + region1.count('C') + \
               region2.count('G') + region2.count('C')
    
    # 返回GC含量（除以2因为每个配对涉及两个碱基）
    return gc_count // 2

def generate_scale_markers():
    """生成刻度标记HTML"""
    return ''.join(
        f'<div style="width: 20px; text-align: center; font-size: 0.75em; color: #6c757d;">'
        f'{i if i % 5 == 0 and i != 0 else ""}</div>'
        for i in range(1, 31)
    )

def generate_base_html(base, annotations, pairings):
    """生成单个碱基的HTML"""
    base_html = f'<div class="base" style="width: 20px; height: 20px; display: flex; ' \
                f'align-items: center; justify-content: center; border: 1px solid #e9ecef; ' \
                f'border-radius: 4px; font-size: 0.9em; font-weight: 500; position: relative;">{base}'
    
    # 添加颜色层
    for ann in annotations:
        hex_color, rgba_color = ann['color']
        layer_id = f"layer-{pairings.index(ann['pairing'])}"
        base_html += f'<div class="color-layer" id="{layer_id}" ' \
                    f'style="position: absolute; top: 0; left: 0; right: 0; bottom: 0; ' \
                    f'background-color: {rgba_color}; mix-blend-mode: multiply; display: none;"></div>'
    
    base_html += '</div>'
    return base_html

def generate_sequence_row(seq, start, end, annotated_chains, chain, pairings):
    """生成一行序列的HTML"""
    row_html = '<div class="sequence-row" style="display: flex; flex-wrap: wrap; margin-bottom: 0px;">'
    for i in range(start, end):
        base = seq[i]
        annotations = [ann for ann in annotated_chains[chain] if ann['start'] <= i <= ann['end']]
        row_html += generate_base_html(base, annotations, pairings)
    row_html += '</div>'
    return row_html

def generate_chain_card(chain, seq, annotated_chains, pairings):
    """生成单个链的卡片HTML"""
    card_html = f"""
    <div class="card" style="margin-bottom: 10px; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1);">
        <div class="card-header" style="background: #f8f9fa; padding: 1px 15px; border-bottom: 1px solid #e9ecef;">
            <span style="font-weight: 600; color: #2c3e50;">{chain}</span>
            <span style="float: right; color: #6c757d; font-weight: 600;">{len(seq)} nt (3'-5') </span>
        </div>
        <div class="card-body" style="padding: 10px;">
            <!-- 序列显示 -->
            <div class="sequence-bases">
                {generate_sequence_row(seq, 0, len(seq), annotated_chains, chain, pairings)}
            </div>
        </div>
    </div>
    """
    return card_html

def read_sequences_from_file(filename):
    """从文件中读取序列信息和比对模式"""
    chains = {}
    mode = 'strict'  # 默认模式
    with open(filename, 'r') as f:
        # 读取第一行判断模式
        first_line = f.readline().strip().lower()
        if first_line in ['strict', 'relaxed']:
            mode = first_line
        else:
            # 如果第一行不是模式指定，则作为序列处理
            if first_line and not first_line.startswith('#'):
                if ':' in first_line:
                    name, seq = first_line.split(':', 1)
                    chains[name.strip()] = seq.strip()
                else:
                    default_name = f"Chain{len(chains)+1}"
                    chains[default_name] = first_line.strip()
        
        # 继续读取剩余行
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                if ':' in line:
                    name, seq = line.split(':', 1)
                    chains[name.strip()] = seq.strip()
                else:
                    default_name = f"Chain{len(chains)+1}"
                    chains[default_name] = line.strip()
    
    return chains, mode

def initialize_aligner(strict=True):
    """初始化比对器"""
    aligner = PairwiseAligner()
    if strict:
        # 严格参数设置
        aligner.mode = 'local'  # 局部比对
        aligner.match_score = 1  # 匹配得分
        aligner.mismatch_score = -100  # 不匹配得分
        aligner.open_gap_score = -100  # 打开空位罚分
        aligner.extend_gap_score = -100  # 扩展空位罚分
    else:
        # 宽松参数设置
        aligner.mode = 'local'  # 局部比对
        aligner.match_score = 1  # 匹配得分
        aligner.mismatch_score = -1  # 不匹配得分
        aligner.open_gap_score = -5  # 打开空位罚分
        aligner.extend_gap_score = -5  # 扩展空位罚分
    return aligner

# 主函数
def main():
    # 从文件读取序列和模式
    chains, mode = read_sequences_from_file('seq_input.txt')
    
    # 检查是否读取到序列
    if not chains:
        print("错误：未读取到任何序列，请检查 seq_input.txt 文件")
        return
    
    # 打印读取到的序列信息
    print("成功读取以下序列：")
    for name, seq in chains.items():
        print(f"{name}: {seq}")
    print(f"比对模式：{mode}")
    
    # 根据模式初始化比对器
    current_aligner = initialize_aligner(strict=(mode == 'strict'))
    
    # 存储配对信息和颜色
    pairings = []
    colors = {}

    # 遍历所有链对，确保每对链只记录一次
    chain_names = list(chains.keys())
    processed_pairs = set()  # 用于存储已经处理过的链对

    for i in range(len(chain_names)):
        for j in range(i + 1, len(chain_names)):  # 从 i+1 开始，避免重复
            chain1 = chain_names[i]
            chain2 = chain_names[j]
            seq1 = chains[chain1]
            seq2 = chains[chain2]

            # 生成链1的反向互补序列
            rc_seq1 = reverse_complement(seq1)

            # 使用 PairwiseAligner 进行比对
            alignments = current_aligner.align(rc_seq1, seq2)
            if alignments:
                align = alignments[0]  # 取最优比对
                # 获取比对长度
                match_length = align.shape[1]  # 比对长度
                if match_length >= 5:
                    # 获取比对区域
                    target_start = align.aligned[0][0][0]  # 链1的起始位置
                    target_end = align.aligned[0][-1][-1]  # 链1的结束位置
                    query_start = align.aligned[1][0][0]  # 链2的起始位置
                    query_end = align.aligned[1][-1][-1]  # 链2的结束位置

                    # 计算链1上的起始和结束位置
                    chain1_start = len(seq1) - target_end
                    chain1_end = len(seq1) - target_start - 1
                    # 链2上的起始和结束位置
                    chain2_start = query_start
                    chain2_end = query_end - 1

                    # 获取比对得分
                    alignment_score = align.score

                    # 获取比对后的序列
                    aligned_seq1 = str(align[0][0])  # 比对后的第一条序列
                    aligned_seq2 = str(align[0][1])  # 比对后的第二条序列
                    
                    # 提取实际配对的碱基对（忽略gap）
                    paired_bases = [(b1, b2) for b1, b2 in zip(aligned_seq1, aligned_seq2) 
                                  if b1 != '-' and b2 != '-']
                    
                    # 计算GC对数量
                    gc_pairs = calculate_gc_content(seq1, seq2, chain1_start, chain1_end, chain2_start, chain2_end)
                    
                    # 修改配对信息记录
                    pairing = (chain1, chain2, chain1_start, chain1_end, chain2_start, chain2_end, 
                               alignment_score, calculate_gc_content(seq1, seq2, chain1_start, chain1_end, chain2_start, chain2_end))
                    pairings.append(pairing)
                    processed_pairs.add(frozenset({chain1, chain2}))

    # 生成颜色
    color_list = generate_colors(len(pairings))
    for i, pairing in enumerate(pairings):
        colors[pairing] = color_list[i]

    # 为每条链创建标注列表
    annotated_chains = {chain: [] for chain in chains}
    for pairing, color in colors.items():
        chain1, chain2, c1_start, c1_end, c2_start, c2_end, score, gc_pairs = pairing
        annotated_chains[chain1].append({'start': c1_start, 'end': c1_end, 'color': color, 'pairing': pairing})
        annotated_chains[chain2].append({'start': c2_start, 'end': c2_end, 'color': color, 'pairing': pairing})

    # 生成 HTML 内容
    html_content = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>IriSeq</title>
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0-alpha1/dist/css/bootstrap.min.css" rel="stylesheet">
        <style>
            body {
                background-color: #ffffff; /* 白色背景 */
                padding: 5px;
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
                color: #1d1d1f; /* 深灰色文字 */
            }
            .container {
                max-width: 1200px;
                margin: 0 auto;
                padding: 5 10px;
            }
            h1 {
                font-size: 1.8rem;
                color: #1d1d1f;
                margin-bottom: 15px;
                text-align: left;  /* 改为左对齐 */
                font-family: 'Trebuchet MS', serif;
                margin-left: 10px;  /* 添加左边距 */
            }
            .card {
                margin-bottom: 10px;
                border: none;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
                background: white;
                transition: transform 0.2s ease, box-shadow 0.2s ease;
            }
            .card:hover {
                transform: translateY(-2px);
                box-shadow: 0 4px 8px rgba(0, 0, 0, 0.15);
            }
            .card-header {
                background: white;
                color: #1d1d1f;
                font-size: 1rem;
                font-weight: 600;
                padding: 5px 15px;
                border-bottom: 1px solid #e0e0e0;
                border-radius: 8px 8px 0 0;
            }
            .sequence-bases {
                display: flex;
                flex-wrap: wrap;
                gap: 3px;
                padding: 10px;
                background-color: white;
                border-radius: 0 0 8px 8px;
            }
            .base {
                padding: 5px 8px;
                border: 1px solid #e0e0e0;
                border-radius: 4px;
                background-color: white;
                font-family: Arial, sans-serif;
                font-weight: bold;
                color: #1d1d1f;
                position: relative;
                font-size: 0.85rem;
            }
            .color-layer {
                position: absolute;
                top: 0;
                left: 0;
                right: 0;
                bottom: 0;
                opacity: 0.6;
                border-radius: 4px;
                display: none; /* 默认隐藏 */
            }
            .controls {
                margin-bottom: 15px;
                padding: 10px;
                background: white;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
            }
            .btn-toggle {
                margin: 3px;
                border-radius: 16px;
                padding: 5px 10px;
                font-size: 0.85rem;
                background-color: white;
                color: #007aff; /* 苹果蓝 */
                border: 1px solid #007aff;
                transition: background-color 0.3s ease, color 0.3s ease;
            }
            .btn-toggle.selected {
                background-color: #007aff;
                color: white;
            }
            .btn-toggle:hover {
                background-color: #007aff;
                color: white;
            }
            .btn-toggle {
                position: relative;
                overflow: hidden;
                transition: all 0.3s ease;
            }
            .score-badge {
                display: inline-block;
                margin-left: 8px;
                padding: 2px 6px;
                border-radius: 12px;
                background-color: rgba(0, 122, 255, 0.1);
                color: #007aff;
                font-size: 0.8em;
                transition: all 0.3s ease;
            }
            .btn-toggle.selected .score-badge {
                background-color: rgba(255, 255, 255, 0.2);
                color: white;
            }
            .btn-toggle:hover .score-badge {
                transform: scale(1.1);
            }
            .color-layer {
                transition: opacity 0.3s ease;
            }
            .gc-badge {
                display: inline-block;
                margin-left: 8px;
                padding: 2px 6px;
                border-radius: 12px;
                background-color: rgba(76, 175, 80, 0.1);
                color: #4CAF50;
                font-size: 0.8em;
                transition: all 0.3s ease;
            }
            .btn-toggle.selected .gc-badge {
                background-color: rgba(255, 255, 255, 0.2);
                color: white;
            }
            .sequence-row {
                display: flex;
                gap: 3px;
                margin-bottom: 5px;
            }
            .base {
                width: 20px;  /* 固定宽度使碱基对齐 */
                height: 20px;
                display: flex;
                align-items: center;
                justify-content: center;
                padding: 0;
            }
            .mode-indicator {
                margin-bottom: 3px;
                padding: 5px 10px;
                border-radius: 16px;
                font-size: 0.9em;
                font-weight: 500;
                display: inline-block;
            }
            .strict-mode {
                background-color: rgba(255, 59, 48, 0.1);
                color: #ff3b30;
            }
            .relaxed-mode {
                background-color: rgba(52, 199, 89, 0.1);
                color: #34c759;
            }
            /* 添加LOGO样式 */
            .header-logo {
                height: 50px;
                margin-right: 50 px;
                margin-bottom: 15px;
                vertical-align: middle;
            }
            .header-container {
                display: flex;
                align-items: center;
                justify-content: flex-start;  /* 改为从左侧开始 */
                padding-left: 10px;  /* 添加左边距 */
            }
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header-container">
                <img src="LOGO.PNG" class="header-logo" alt="IriSeq Logo">
                <h1>
                    <span style="font-weight: bold;">IriSeq</span>: 
                    <span style="font-weight: normal;">Colors reveal pairing</span>
                </h1>
            </div>
            <div class="controls card">
                <div class="card-body">
                    <div class="d-flex flex-wrap gap-2">
    """

    # 根据当前比对器设置模式指示器
    is_strict_mode = current_aligner.mismatch_score == -100
    strict_style = ""  # 移除内联样式，由JavaScript控制
    relaxed_style = ""  # 移除内联样式，由JavaScript控制

    html_content += f"""
            <div class="mode-indicator strict-mode" data-initial-mode="{str(is_strict_mode).lower()}">Strict Mode</div>
            <div class="mode-indicator relaxed-mode">Relaxed Mode</div>
    """

    # 添加带有得分信息的按钮
    for i, pairing in enumerate(pairings):
        hex_color, rgba_color = colors[pairing]
        chain1, chain2, _, _, _, _, score, gc_pairs = pairing
        html_content += f"""
        <button class="btn btn-toggle" onclick="toggleLayer('layer-{i}')" data-score="{score}">
            <span class="pair-names">{chain1} & {chain2}</span>
            <span class="score-badge">{int(score)} bp</span>
            <span class="gc-badge">{gc_pairs} GC pairs</span>
        </button>
        """

    html_content += """
                    </div>
                </div>
            </div>
    """

    # 生成每个链的卡片
    for chain, seq in chains.items():
        html_content += generate_chain_card(chain, seq, annotated_chains, pairings)
    
    # 添加 JavaScript 和结束标签
    html_content += """
            <script>
                // 初始化模式显示
                function initializeMode() {
                    const strictIndicator = document.querySelector('.strict-mode');
                    const relaxedIndicator = document.querySelector('.relaxed-mode');
                    const initialMode = strictIndicator.getAttribute('data-initial-mode') === 'true';
                    
                    if (initialMode) {
                        strictIndicator.style.display = 'inline-block';
                        relaxedIndicator.style.display = 'none';
                    } else {
                        strictIndicator.style.display = 'none';
                        relaxedIndicator.style.display = 'inline-block';
                    }
                }

                // 页面加载时初始化
                window.onload = function() {
                    initializeMode();
                    initializeButtons();
                };

                // 改进的切换效果
                function toggleLayer(layerId) {
                    const layers = document.querySelectorAll(`.color-layer#${layerId}`);
                    const button = document.querySelector(`button[onclick*="${layerId}"]`);
                    
                    // 切换按钮状态
                    button.classList.toggle('selected');
                    
                    // 添加更流畅的过渡效果
                    layers.forEach(layer => {
                        if (layer.style.display === 'none') {
                            layer.style.display = 'block';
                            layer.style.opacity = 0;
                            setTimeout(() => {
                                layer.style.opacity = 1;
                            }, 10);
                        } else {
                            layer.style.opacity = 0;
                            setTimeout(() => {
                                layer.style.display = 'none';
                            }, 300);
                        }
                    });
                }

                // 改进的按钮初始化
                function initializeButtons() {
                    const buttons = document.querySelectorAll('.btn-toggle');
                    buttons.forEach(button => {
                        const layerId = button.getAttribute('onclick').match(/'(.*?)'/)[1];
                        const layers = document.querySelectorAll(`.color-layer#${layerId}`);
                        const isVisible = layers[0].style.display !== 'none';
                        button.classList.toggle('selected', isVisible);
                    });
                }

                // 添加按钮点击的涟漪效果
                document.addEventListener('click', function(e) {
                    if (e.target.classList.contains('btn-toggle')) {
                        const rect = e.target.getBoundingClientRect();
                        const x = e.clientX - rect.left;
                        const y = e.clientY - rect.top;
                        
                        const ripple = document.createElement('span');
                        ripple.className = 'ripple-effect';
                        ripple.style.left = `${x}px`;
                        ripple.style.top = `${y}px`;
                        e.target.appendChild(ripple);
                        
                        setTimeout(() => {
                            ripple.remove();
                        }, 600);
                    }
                });
            </script>
        </div>
        <!-- 添加页尾 -->
        <footer style="
            margin-top: 20px;
            padding: 20px;
            text-align: center;
            color: #6c757d;
            font-size: 0.9em;
            border-top: 1px solid #e9ecef;
        ">
            for the aid of understanding, design and screening
        </footer>
    </body>
    </html>
    """

    # 将 HTML 写入文件
    with open('dna_alignment_visualization.html', 'w') as f:
        f.write(html_content)

    print("HTML 文件已生成：dna_alignment_visualization.html")

if __name__ == "__main__":
    main()