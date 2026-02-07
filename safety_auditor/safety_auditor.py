def audit_safety(context, report, ask_llm):

    try:
        # STEP 6: SAFETY AUDIT
        check_prompt = f"""
        You are a scientific validation auditor.

        ### TASK:
        Compare the REPORT against the SOURCE DATA.

        ### SOURCE DATA:
        <context>{context}</context>

        ### REPORT TO CHECK:
        {report}

        ### RULES:
        - Every factual claim in the REPORT must be explicitly supported by the SOURCE DATA.
        - Claims must match the described mechanisms, results, and scope.
        - If a claim is missing, inferred, exaggerated, or not directly stated in the SOURCE DATA -> OVERCLAIMED.
        - Minor rewording is allowed; new meaning is not.
        - If ALL claims are supported -> SAFE.

        ### OUTPUT:
        Return ONLY one word:
        SAFE or OVERCLAIMED
        """
        safety = ask_llm(check_prompt)

        # Display Final Report
        print("\n" + "="*60)
        print("ANALYSIS REPORT")
        print("="*60 + "\n")
        print(report)
        print("\n" + "="*60)
        
        # Warning
        if "OVERCLAIMED" in safety:
            print("\n[SAFETY NOTICE]: Some claims may exceed the available evidence. Please verify with primary sources.")
                
    except Exception as e:
        print(f"\n[ERROR] {str(e)}")

            