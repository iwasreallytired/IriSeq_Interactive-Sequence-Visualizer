# IriSeq (Interactive Sequence Visualizer) 
IriSeq is a visualization tool specifically designed for nucleic acid sequence analysis. Its goal is to assist researchers in intuitively observing and analyzing pairing interactions between nucleic acid sequences through efficient algorithms and an interactive interface. With a visually appealing interface and precise alignment capabilities as its core features, the software integrates functions such as pairing analysis, GC content calculation, and color annotation, making it ideal for visualizing base complementarity results in nucleic acid research.

## Principle: 
The alignment between nucleic acid sequences is implemented using BioPython's PairwiseAligner. The local alignment algorithm (Smith-Waterman algorithm) is employed to extract the most similar local regions from the sequences，which is
- set to local, focusing only on similar regions within the sequences.
- determined based on the scoring settings for matches, mismatches, gap opening, and gap extension.

## Features:
1. Dual Alignment Modes
- Strict Mode: Employs high penalty scores, suitable for highly specific complementary pairings.
- Relaxed Mode: Uses low penalty scores, ideal for screening potential non-orthogonal sequences.
2. Color Generation - Utilizes HTML and CSS technologies to create an interactive interface for displaying sequence alignment results. Dynamic style adjustments (e.g., layered color effects) enhance the user experience.
3. GC Content Calculation - Provides GC content information to assist users in analyzing sequence stability.

## Scenarios:
1. Teaching Demonstrations - Analyze the similarity and predict the basic process of sequence pairing.
2. Experimental Reproduction - Visualize paired regions to aid in understanding reported DNA hybridization reactions.
3. Sequence Design - Assist in designing DNA sequences to ensure specific pairings between sequences.
4. Library Screening: Facilitate the evaluation of the orthogonality and stability of large-scale DNA libraries.
### Additional Benefits of Combining IriSeq with NUPACK:
  - Comprehensive Sequence Analysis - Integrating pairing region information with secondary structure predictions for a deeper understanding of sequence characteristics.
  - Enhanced Design Efficiency - IriSeq offers rapid screening, while NUPACK provides detailed verification. This combination is particularly suited for DNA nanostructure design, multi-strand interaction analysis, and hybridization chain reaction analysis.



