# Design gRNA

CRISPR/Cas9用のgRNAを見繕うパイプライン

1. 標的領域を入力に、gRNA候補配列を列挙する（PAM配列は、`NGG`にのみ対応）
2. <a href="https://github.com/snugel/cas-offinder">Cas-OFFinder</a> (Bae S, _et al_., _Bioinformatics_, 2014)を用いて、各gRNA候補によるオフターゲット可能性がある配列をリファレンス配列から検索する
3. <a href="https://flycrispr.org/target-finder">flyCRISPR</a>のOptimal Target Finder (OTF)の評価基準に準拠して、gRNA候補配列を絞る

## 経緯

OTFの評価基準に準じているなら、OTFを使えばいいじゃない、と思われるかもしれないので、このパイプラインを作成した経緯を記す。<br>
CRISPRiを用いて、D. melのレトロトランスポゾン（TEs）の発現を低下させる実験をしたかった。しかし、ショウジョウバエでは意外とCRISPRiの技術、特に発現抑制に関する技術が少なかった。調べてみると、tgCRISPRi <a href="https://www.nature.com/articles/s41467-023-40836-3">(Auradkar, _et al_., _Nat. Commun_., 2023)</a>という技術が存在していたので、これを採用しようと思った。<br>そこで、OTFを使って、tgRNAをデザインしようとしたが、2つの問題が生じた。1つは、OTFでは16 nt以下のgRNA長を指定できないこと、もう1つは、OTFではショウジョウバエゲノム（dm6）をリファレンスにしているが、TEsはゲノム中に複数コピーあるので、必然的にどのgRNA候補もオフターゲットの可能性が高いと表示されたこと。<br>
これらを解決するために、このパイプラインを作成した。

したがって、このパイプラインはOTFではできない次のことができる。

- 任意のリファレンス配列を指定できる（OTFで指定できるのはD. melを中心としたショウジョウバエゲノムのみ）
- gRNA（プロトスペーサー配列）長を短く指定できる（OTFで指定できるのは16-20 ntのみ）

なお、リファレンスゲノムの問題については、`dm6_singlecopyTEs.fa` <a href="https://genome.cshlp.org/content/suppl/2020/10/14/gr.259200.119.DC1">(Treiber & Waddell, *Genome Res*, 2020)のsupplement</a>を使って解決した。<br>
これは、dm6をベースに、ゲノム中の反復配列にマスクをかけ、各TEsを別chromosome扱いで追記したfastaファイル。<br>
レアかもしれないが、同じ境遇で困った人がいたら、参考にしてください。

## OTFのオフターゲット評価に準じた絞り込み

gRNA配列を次の2領域に分け、

- Proximal seq: PAM配列近傍の12 nt
- Distal seq: 残りの配列

（例）GATCGATGTTGAATGCCGATTGG (20 nt gRNA)とあったら、
|Distal|Proximal|PAM|
|--|--|--|
|GATCGATG|TTGAATGCCGAT|TGG|

次の条件に当てはまるオフターゲットサイトを1つでも持つgRNA候補を除外する。

```
high stringencyの場合：

- Proximal seqに0 mismatch
- Proximal seqに1 mismatch かつ Distal seqに2 mismatch未満

のいずれかに該当する
```

```
maximum stringencyの場合：

- Proximal seqに0 mismatch
- Proximal seqに1 mismatch かつ Distal seqに5 mismatch未満
- Proximal seqに2 mismatch かつ Distal seqに2 mismatch未満

のいずれかに該当する
```

ただし、オンターゲットが全長0ミスマッチとなるので、全長0ミスマッチサイトについてはすべて許容するようにしている。<br>
全長0ミスマッチサイトには、出力するテーブルの`mismatch`カラムに`# 0MM`フラグを立てている。<br>
gRNA長が短いと、稀に全長0ミスマッチでもオフターゲットであることがあるので、「`# 0MM`フラグが立っているがオンターゲットでない」サイトが含まれていないかを必ず確認する。

## 使用例

```
python3 gRNA_finder.py -i input_file.txt -g reference_file.fa
```

## 入出力

- 入力は、プレーンテキスト / fasta
- 出力は、tsv

## オプション

|オプション名|内容|
|--|--|
|`-i`, `--input_file`|入力ファイルのパスを指定する、必須|
|`-l`, `--gRNA_length`|PAM配列より前の配列の長さを指定する、数値、デフォルト値=`20`|
|`-g`, `--genome_path`|リファレンスファイルのパス、fasta / 2BITファイルに対応（Cas-OFFinderの要件）<br>必須|
|`-m`, `--mismatches`|ミスマッチ数を指定する、数値、デフォルト値=`3`|
|`-b`, `--bulge`|バルジサイズを指定する、数値、デフォルト値=`0`|
|`--CO_path`|Cas-OFFinderのパスを指定、デフォルト値=`cas-offinder`|
|`-s`, `--stringency`|stringencyのレベルを`high`または`maximum`で指定、デフォルト値=`high`|
|`-t`, `--save_temp`|Cas-OFFinderの入力ファイルおよび出力ファイルを保存するかどうか<br>（このオプションをつけない場合は保存しない）|
|`-o`, `--output_file`|出力ファイルのパスを指定する、デフォルト値は`output_GF.txt`|

## 必要な要件

- python3
- Cas-OFFinder (3.0.0を想定)

が使えること。
