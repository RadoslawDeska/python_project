def compare_solvents(current_solvent, solvent_memory):
    return current_solvent.lower() in [memo.lower() for memo in solvent_memory]