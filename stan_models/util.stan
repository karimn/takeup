int num_test(array[] int to_test, array[] int target_val, int test_equality) {
  int num_to_test = num_elements(to_test);
  int num_targets = num_elements(target_val);
  int result = 0;
  
  array[num_to_test] int sorted_to_test = sort_asc(to_test);
  
  for (to_test_index in 1:num_to_test) {
    int found = 0;
    
    for (target_index in 1:num_targets) {
      if (sorted_to_test[to_test_index] == target_val[target_index]) {
        if (test_equality) {
          result += 1;
        }
        
        found = 1;
        break;
      } 
    }
    
    if (!found && (1 - test_equality)) {
      result += 1;
    }
  }
  
  return(result);
}

int num_equals(array[] int to_test, array[] int target_val) {
  return(num_test(to_test, target_val, 1));
}

array[] int count(int count_size, array[] int find_in) {
  array[count_size] int count_array = rep_array(0, count_size);
  
  for (count_index in 1:count_size) {
    count_array[count_index] = num_equals(find_in, { count_index });
  }
  
  return(count_array);
}

array[] int count_by_group_test(array[] int to_count, array[] int group, array[] int target_val, int test_equality) {
  int num_to_count = num_elements(to_count);
  int num_groups = max(group); // num_elements(unique(group));
  array[num_groups] int group_sizes = count(num_groups, group); 
  array[num_to_count] int to_count_group_sorted = to_count[sort_indices_asc(group)];
  
  array[num_groups] int group_count = rep_array(0, num_groups);
  int group_pos = 1;
  
  if (num_elements(group) != num_to_count) {
    reject("Incompatible array sizes.");
  }
  
  for (group_index in 1:num_groups) {
    if (group_sizes[group_index] > 0) {
      int group_end = group_pos + group_sizes[group_index] - 1;
      
      group_count[group_index] = num_test(to_count[group_pos:group_end], target_val, test_equality);
      
      group_pos = group_end + 1;
    }
  }
  
  return(group_count);
}

array[] int which(array[] int to_test, array[] int target_val, int test_equality) {
  int num_to_test = num_elements(to_test);
  int num_targets = num_elements(target_val);
  int num_which = num_test(to_test, target_val, test_equality);
  array[num_which] int result;
  int curr_result_index = 0;
  
  for (to_test_index in 1:num_to_test) {
    int found = 0;
    
    for (target_index in 1:num_targets) {
      if (to_test[to_test_index] == target_val[target_index]) {
        if (test_equality) {
          curr_result_index += 1;
          result[curr_result_index] = to_test_index;
        }
       
        found = 1; 
        break;
      } 
    }
        
    if (!found && (1 - test_equality)) {
      curr_result_index += 1;
      result[curr_result_index] = to_test_index;
    }
    
    if (curr_result_index >= num_which) {
      break;
    }
  }
  
  return(result);
}

array[] int seq(int from, int to, int by) {
  int reverse = from > to;
  int seq_len = (((1 - 2 * reverse) * (to - from)) %/% by) + 1;
  array[seq_len] int result_seq;
  
  for (seq_index in 1:seq_len) {
    result_seq[seq_index] =  from + (1 - 2 * reverse) * ((seq_index - 1) * by);
  }
  
  return(result_seq);
}

array[] int array_subtract(array[] int left, array[] int right) {
  int array_size = num_elements(left);
  array[array_size] int array_sub;
  int right_array_size = num_elements(right);  
  
  if (right_array_size != array_size && right_array_size != 1) {
    reject("Incompatible array sizes.");
  }
  
  for (array_index in 1:array_size) {
    array_sub[array_index] = left[array_index] - right[right_array_size > 1 ? array_index : 1];
  }
  
  return(array_sub);
}

array[] int array_product(array[] int left, array[] int right) {
  int array_size = num_elements(left);
  array[array_size] int array_prod;
  int right_array_size = num_elements(right);  
  
  if (right_array_size != array_size && right_array_size != 1) {
    reject("Incompatible array sizes.");
  }
  
  for (array_index in 1:array_size) {
    array_prod[array_index] = left[array_index] * right[right_array_size > 1 ? array_index : 1];
  }
  
  return(array_prod);
}

array[] int array_add(array[] int left, array[] int right) {
  int array_size = num_elements(left);
  array[array_size] int array_sum;
  int right_array_size = num_elements(right);  
  
  if (right_array_size != array_size && right_array_size != 1) {
    reject("Incompatible array sizes.");
  }
  
  for (array_index in 1:array_size) {
    array_sum[array_index] = left[array_index] + right[right_array_size > 1 ? array_index : 1];
  }
  
  return(array_sum);
}

int in_array(int to_test, array[] int target_val) {
  int num_target = num_elements(target_val);
  
  for (target_index in 1:num_target) {
    if (to_test == target_val[target_index]) return(1);
  }
  
  return(0);
}