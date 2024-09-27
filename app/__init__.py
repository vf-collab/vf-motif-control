# app/__init__.py

from .codon_optimizer import codon_optimization
from .sequence_utils import not_in, tuple_type, validate_input, GC, check_correct, find_expansion, opt_expansion, repeat_checker, find_hairpin_repeat, prime, back_translate, repeat_remove, hairpin_remove, screen_sequence, pattern_generator
from .enzyme_utils import is_CpG_island, motif_remove, CpG_check, CpG_remove, CpG_island_remove
from .progress_utils import update_progress_bar
